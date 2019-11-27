#!/usr/bin/env julia

## test_network_lending.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test NetworkLending and plot frequency trajectories.

using Revise
using NetworkLending
using PyPlot

N = 1000 # population size
k = 100 # network degree
r = 1.5 # synergy factor
z = 1.0 # interest
d = 0.9 # decrement for bad actors
w = 0.01 # selection strength

good_reputations = [1, 3]

# initialize the game and network basics
gp = GameParams(r, z, d)
game = NetworkGame(w, gp, good_reputations)
network = Network(N, k)

# initialize the NetworkPopulation object
pop = NetworkPopulation(network, game)

# initialize the frequency trajectories
ft = FreqTraj(pop)

for i in 1:N
	evolve_and_track!(pop, ft)
end

function strategy_label(
	strategy::Int64
	)
	# return the label for an individual's strategy
	strat_string = ""
	strat_string *= (strategy ÷ 2 == 0 ? "defect" : "cooperate")
	strat_string *= ", "
	strat_string *= (strategy % 2 == 0 ? "shirk" : "payback")
	return strat_string
end

fig = plt.figure()
[plt.plot(ft.freqs[:,x], label=("$(strategy_label(x-1))")) for x in 1:4]
legend(loc=2)
plt.xlabel("time")
plt.ylabel(L"\rho")
plt.title("all frequencies")
plt.tight_layout()
display(fig)

#
# fig = plt.figure()
# [plt.plot(payback_freqs[x,:], label=("$(x-1)")) for x in 1:2]
# legend(loc=2)
# plt.title("payback")
# plt.tight_layout()
# display(fig)
#
#
# fig = plt.figure()
# [plt.plot(coop_freqs[x,:], label=("$(x-1)")) for x in 1:2]
# legend(loc=2)
# plt.title("coop")
# plt.tight_layout()
# display(fig)
