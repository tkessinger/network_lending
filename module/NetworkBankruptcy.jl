#!/usr/bin/env julia

## NetworkBankruptcy.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Module for simulating central bank lending and repayment
## on a regular network.
## This differs from NetworkLending in the following ways:
## 1. There are only two strategies, cooperate and defect.
## 2. Banks punish individuals whose fitness falls below some threshold
## by lending them less money in subsequent rounds.
## For now the threshold is zero.

module NetworkBankruptcy

    using Random, GenRRG, StatsBase

    export Network, NetworkGame, GameParams, NetworkPopulation
	export FreqStats, FreqTraj
    export evolve!, evolve_and_track!, evolve_and_track_until_fixation!
	export pair_PGG_payoffs

	struct Network
		# static type for storing the graph itself

		N::Int64 # number of nodes
		k::Int64 # degree of network
		neighbors::Array{Array{Int64, 1}, 1} # each individual's neighbors
		edges::Array{Array{Int64, 1}, 1} # all edges on the network
		adjacency::Array{Bool, 2} # adjacency matrix, true if connected

		# constructors
		function Network(N::Int64, k::Int64)
			# this uses GenRRG's graph generation function
			# to create a connected RRG from scratch
			neighbors, edges, adjacency = generate_connected_RRG(N, k)
			return new(N, k, neighbors, edges, adjacency)
		end
		function Network(N::Int64, k::Int64, neighbors::Array{Array{Int64, 1}}, edges::Array{Array{Int64, 1}}, adjacency::Array{Bool, 2})
			# this is just in case you already have a graph in mind
			return new(N, k, neighbors, edges, adjacency)
		end
	end

    struct GameParams
        # new type for "game" parameters
        # this should minimize the amount of passing around that needs to be done
		# the game is a pairwise PGG with lending

        r::Float64 # synergy parameter
        d::Float64 # decrement for individuals with bad reputations

        # constructors
        function GameParams(r::Float64)
			# initializes a game with no decrement
            return new(r, 0.0)
        end
        function GameParams(r::Float64, d::Float64)
            return new(r, d)
        end
    end

	struct NetworkGame
		# type for storing game and miscellaneous parameters

        w::Float64 # strength of selection
        params::GameParams # parameters of the PGG with lending
		thresholds::Array{Float64,1} # the payback thresholds for each game
		# set this to [1.0-d, 1.0] to ensure that people must pay back the principal
		function NetworkGame(w::Float64, params::GameParams)
			# initializes a game with default reputations
			return new(w, params, [1.0, 1.0-params.d])
		end
		function NetworkGame(w::Float64, params::GameParams, thresholds::Array{Float64, 1})
			return new(w, params, thresholds)
		end
    end

    mutable struct NetworkPopulation
        # the population of individuals and all information about them
        # this should make it easier to pass around neighbors, etc.

		network::Network # not mutable: the network
        game::NetworkGame # not mutable: game attributes
		strategies::Array{Int64, 1} # array of strategies
		reputations::Array{Int64, 1} # array of reputations
		fitnesses::Array{Float64, 1} # array of fitnesses
        generation::Int64
        verbose::Bool # turn this on for error tracking

        # constructor if game and network are already specified
        function NetworkPopulation(
            network::Network,
			game::NetworkGame
            )
			# begin by initializing the population with random strategies
			strategies = random_strategies(network.N, [0,1])
			# update everyone's reputations and fitnesses
			#reputations = zeros(Int64, network.N)
			reputations = [rand([0,1]) for i in 1:network.N]
			fitnesses = get_all_fitnesses(network, game, strategies, reputations)
            generation = 0
			return new(network, game, strategies, reputations, fitnesses, generation)
        end
    end

	struct FreqStats
		# type for storing instantaneous frequency statistics

		freqs::Array{Float64, 1} # frequency of each strategy x_i
		coop_freqs::Array{Float64, 1} # cooperator frequency
		payback_freqs::Array{Float64, 1} # payback frequency
		pair_freqs::Array{Float64, 2} # pair frequencies x_{ij}
		conditional_freqs::Array{Float64, 2} # conditional frequencies q_{i|j}

		# constructors
		function FreqStats(
			# constructor in case frequencies have already been obtained
			freqs::Array{Int64, 1},
			coop_freqs::Array{Int64, 1},
			payback_freqs::Array{Int64, 1},
			pair_freqs::Array{Int64, 2},
			conditional_freqs::Array{Int64, 2}
			)
			return new(freqs, coop_freqs, payback_freqs, pair_freqs, conditional_freqs)
		end
		function FreqStats(pop::NetworkPopulation)
			# constructor if frequencies haven't been obtained
			# this obtains them first
			freqs = zeros(Float64, 4)
			coop_freqs = zeros(Float64, 2)
			payback_freqs = zeros(Float64, 2)
			pair_freqs = zeros(Float64, 4, 4)
			conditional_freqs = zeros(Float64, 4, 4)

			for indv in 1:pop.network.N
				freqs[pop.strategies[indv]+1] += 1.0
			end
			freqs /= pop.network.N
			# types 0 and 2 are defectors, 1 and 3 are cooperators
			coop_freqs[1] = freqs[1] + freqs[3]
			coop_freqs[2] = freqs[2] + freqs[4]
			# types 0 and 1 are shirkers, 2 and 3 are payersback
			payback_freqs[1] = freqs[1] + freqs[2]
			payback_freqs[2] = freqs[3] + freqs[4]

			# pair and conditional frequencies
			for indv in 1:pop.network.N
                for neighbor in pop.network.neighbors[indv]
                    pair_freqs[pop.strategies[indv]+1, pop.strategies[neighbor]+1] += 1.0
                end
            end
            # pair frequencies should be symmetric: x_{ij} = x_{ji}
            pair_freqs += transpose(pair_freqs)
            # there are Nk/2 total edges
            pair_freqs /= (2*pop.network.N*pop.network.k)
			# q_{i|j} = x_{ij}/x_j
			conditional_freqs = pair_freqs./transpose(freqs)

			return new(freqs, coop_freqs, payback_freqs, pair_freqs, conditional_freqs)
		end
	end

	mutable struct FreqTraj
		# type for storing the time trajectories of each strategy, etc.

		freqs::Array{Float64, 2} # frequency of each strategy x_i
		coop_freqs::Array{Float64, 2} # cooperator frequency
		payback_freqs::Array{Float64, 2} # payback frequency
		pair_freqs::Array{Float64, 3} # pair frequencies x_{ij}
		conditional_freqs::Array{Float64, 3} # conditional frequencies q_{i|j}

		# constructor
		function FreqTraj(pop::NetworkPopulation)
			# the first index of each array is the current generation
			# these are moran generations
			# basic idea is that at each time step, a new FreqStats instance is ``appended''

			freqs = zeros(Float64, 1, 4)
			coop_freqs = zeros(Float64, 1, 2)
			payback_freqs = zeros(Float64, 1, 2)
			pair_freqs = zeros(Float64, 1, 4, 4)
			conditional_freqs = zeros(Float64, 1, 4, 4)

			fs = FreqStats(pop)
			freqs[1,:] = fs.freqs
			coop_freqs[1,:] = fs.coop_freqs
			payback_freqs[1,:] = fs.payback_freqs
			pair_freqs[1,:,:] = fs.pair_freqs
			conditional_freqs[1,:,:] = fs.conditional_freqs

			return new(freqs, coop_freqs, payback_freqs, pair_freqs, conditional_freqs)
		end
	end

	function ugly_append(
		array1::Array{Float64, 3},
		array2::Array{Float64, 2}
		)
		# hack for appending a 2d array to a 3d array
		# for FreqTraj time series
		# weirdly, cat() does not seem to allow this
		tmp_size = size(array1)
		tmp_array = zeros(Float64, tmp_size[1]+1, tmp_size[2], tmp_size[3])
		tmp_array[1:end-1, :, :] = array1
		tmp_array[end, :, :] = array2
		array1 = tmp_array
		return array1
	end

	function track_freqs!(
		pop::NetworkPopulation,
		ft::FreqTraj
		)
		# updates a FreqTraj object with the current FreqStats

		fs = FreqStats(pop)

		ft.freqs = vcat(ft.freqs, fs.freqs')
		ft.coop_freqs = vcat(ft.coop_freqs, fs.coop_freqs')
		ft.payback_freqs = vcat(ft.payback_freqs, fs.payback_freqs')

		# see above hack for appending a 2d array to a 3d array
		ft.pair_freqs = ugly_append(ft.pair_freqs, fs.pair_freqs)
		ft.conditional_freqs = ugly_append(ft.conditional_freqs, fs.conditional_freqs)
	end

	function random_strategies(
		n::Int64,
		strategy_set::Array{Int64,1}=collect(0:1)
		)
		# initializes a population with random strategies
		# if n strategies are assigned,
		# then N/n individuals will have each strategy
		strategies = zeros(Int64, n)
		permutation = randperm(n)
		g = length(strategy_set)
		for (si, strat) in enumerate(strategy_set)
			[strategies[x] = strat for x in (si-1)*n÷g+1:si*n÷g]
		end
		return strategies
	end

    function get_all_fitnesses(
		network::Network,
		game::NetworkGame,
		strategies::Array{Int64, 1},
        reputations::Array{Int64, 1},
        verbose::Bool=false
        )
		# gets the fitnesses of the entire population
		# first computes payoffs, then uses w to obtain fitnesses
        payoffs = zeros(Float64, network.N)
		fitnesses = ones(Float64, network.N)
		# iterating over pairs saves some time here
		# and avoids double counting
        for (pi, pair) in enumerate(network.edges)
            payoffs[pair] += pair_PGG_payoffs(strategies[pair], reputations[pair], game.params, verbose)
        end
		# f_i = 1 - w + w*p_i
		[fitnesses[x] += game.w*(payoffs[x] - 1.0) for x in 1:network.N]
        return fitnesses
    end

	function get_all_fitnesses(
		pop::NetworkPopulation
		)
		# a wrapper for the above, with a NetworkPopulation object
		return get_all_fitnesses(pop.network, pop.game, pop.strategies, pop.reputations, pop.verbose)
	end

	function obtain_indv_payoff(
		pop::NetworkPopulation,
		indv::Int64
		)
		# returns the payoff for a single individual
		# based on their neighbors
		payoff = 0
		for (ni, neighb) in enumerate(pop.network.neighbors[indv])
			pair = [indv, neighb]
			payoff += pair_PGG_payoffs(pop.strategies[pair], pop.reputations[pair], pop.game.params, pop.verbose)[1]
		end
		return payoff
	end

	function update_indv!(
		pop::NetworkPopulation,
		invader::Int64,
		invadee::Int64
		)
		# updates both the strategy and fitness of an individual
		# note that update_reputations!() is NOT called here
		# we may want to keep that separate
		update_strategies!(pop, invader, invadee)
		# we need to update both the invadee's fitness
		# and that of all their neighbors
		update_fitness!(pop, invadee)
		for (ni, neighb) in enumerate(pop.network.neighbors[invadee])
			update_fitness!(pop, neighb)
		end
	end

	function update_strategies!(
		pop::NetworkPopulation,
		invader::Int64,
		invadee::Int64
		)
		# changes an individual's strategy
		# to that of an invading neighbor
		pop.strategies[invadee] = pop.strategies[invader]
	end

	function update_fitness!(
		pop::NetworkPopulation,
		indv::Int64
		)
		# updates a single individual's payoff
		payoff = obtain_indv_payoff(pop, indv)
		pop.fitnesses[indv] = 1.0 - pop.game.w + pop.game.w*payoff
	end

    function update_all_fitnesses!(
        pop::NetworkPopulation
        )
		# updates the entire population's payoffs
        pop.fitnesses = get_all_fitnesses(pop)
    end

	function update_all_reputations!(
        pop::NetworkPopulation
        )
		# updates the entire population's reputations
        for indv in 1:pop.network.N
            update_reputation!(pop, indv)
		end
    end

	function update_reputation!(
		pop::NetworkPopulation,
		indv::Int64
		)
		# updates a single individual's reputation
		# depending on their current strategy
		# this is just a wrapper for get_reputation() below
		pop.reputations[indv] = get_reputation(pop.fitnesses[indv], pop.reputations[indv], pop.game.thresholds, pop.network.k)
	end

	function get_all_reputations(
		pop::NetworkPopulation
		)
		return [get_reputation(pop.fitnesses[i], pop.reputations[i], pop.game.thresholds, pop.network.k) for i in 1:pop.network.N]
	end

	function get_reputation(
		fitness::Float64,
		reputation::Int64,
		thresholds::Array{Float64,1},
		k::Int64
		)
		# if the individual has enough "money" to pay back their k loans, return 1
		# else, return 0
		println("indv has fitness $fitness and reputation $reputation")
		println("this is $(fitness > k*thresholds[reputation+1] ? "higher" : "lower") than the threshold $(k*thresholds[reputation+1])")
		println("so individual's reputation is set to $(fitness > k*thresholds[reputation+1] ? 1 : 0)")
		if fitness > k*thresholds[reputation+1]
			return 1
		else
			return 0
		end
	end

    function pair_PGG_payoffs(
        strat::Array{Int64, 1},
        rep::Array{Int64, 1},
        gp::GameParams,
        verbose::Bool=false
        )
		# computes the pair payoffs for a PGG with lending
		# each individual is given 1 unit of capital
		# individuals with bad reputations are given (1-d) instead
		# they invest or do not invest it in a common pairwise pool
		# pool amount is multiplied by r and redistributed
		# the bank demands z back as interest

        r, d = gp.r, gp.d # get the game parameters

		# do these individuals have good reputations?
        money = [(x==1 ? 1.0 : 1.0-d) for x in rep]

		# get the individual strategy bits
		coop = strat

		# figure out how much money goes in the pot
        PGG_pot = r*(sum([coop[i]*money[i] for i in 1:2]))

		# determine payoffs
		payoffs = [PGG_pot/2 + (-coop[i])*money[i] for i in 1:2]

        return payoffs
    end

    function evolve!(
		pop::NetworkPopulation,
		num_gens::Int64=1,
		verbose::Bool=false
		)
		# evolve the population for a specified number of generations
		# using death-birth updating
		for i in 1:num_gens
			# choose a random individual to die
			indv = rand(1:pop.network.N)
			# choose a random neighbor to replace them
			# weighted by fitnesses
			neighbor_fitnesses = pop.fitnesses[pop.network.neighbors[indv]]
			invading_neighbor = sample(pop.network.neighbors[indv], Weights(neighbor_fitnesses))
			# if any strategies have changed, update strategy and fitnesses
			if pop.strategies[invading_neighbor] != pop.strategies[indv]
				update_indv!(pop, invading_neighbor, indv)
			end
			# update their reputation
			# later we can allow neighbors to update their reputations, too
			update_reputation!(pop, indv)
			pop.generation += 1
		end
	end

	function evolve_and_track!(
		pop::NetworkPopulation,
		ft::FreqTraj
		)
		# evolves a population N time steps (one moran generation)
		# and updates the frequency trajectories
		evolve!(pop, pop.network.N)
		track_freqs!(pop, ft)
	end

	function evolve_and_track_until_fixation!(
		pop::NetworkPopulation,
		ft::FreqTraj
		)
		while !any([freq == 1.0 for freq in ft.freqs[end,:]])
			evolve_and_track!(pop, ft)
		end
	end

# final end statement to close the module
end
