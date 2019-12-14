#!/usr/bin/env julia

## parallel_lending.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Parallelized implementation of NetworkLending.
## Obtain mean type frequencies for different values of r and z.

using Random, Statistics
using NetworkLending
using Distributed
using Revise
using ArgParse
using CSV
using Dates
using DataFrames
import JSON

function read_parameters(defpars::Dict{String,Any},
    inputfile = nothing)

    pars = copy(defpars)

    # read JSON file
    if inputfile != nothing
        inpars = JSON.parsefile(inputfile)
    else
        inpars = Dict()
    end

    for parkey in keys(defpars)
        if "type" in keys(pars[parkey])
            if isprimitivetype(pars[parkey]["type"]) ||
                pars[parkey]["type"] == String
                T = pars[parkey]["type"]
            end
        else
            # default type is Float64
            T = Float64
        end
        #println(parkey, T)
        if T <: Int
            convertf = (val)->round(T, val)
        else
            convertf = (val)->convert(T, val)
        end

        # use defpars for list of usable parameters in JSON
        if parkey in keys(inpars)
            if "value" in keys(inpars[parkey])
                val = inpars[parkey]["value"]
            elseif "range" in keys(inpars[parkey])
                valr = inpars[parkey]["range"]
                if "log" in keys(valr)
                    b = valr["log"]
                    rf = (r)->b.^r
                    pop!(valr, "log")
                else
                    rf = (r)->r
                end
                start = pop!(valr, "start")
                rkws = Dict(zip(Symbol.(keys(valr)), values(valr)))
                val = rf(range(start; rkws...))
            end
        else
            val = pars[parkey]["value"]
        end

        if !isstructtype(typeof(val)) || typeof(val) == String || typeof(val) == Bool
            pars[parkey] = [convertf(val)]
        else
            pars[parkey] = convertf.(val)
        end
    end

    return pars
end

function main(args)

    s = ArgParseSettings(description =
        "run NetworkLending simulations across multiple cores")
    @add_arg_table s begin
        "--ncpus"
            arg_type = Int64
            default = max(round(Int, Sys.CPU_THREADS / 2), 1)
        "--input"
            default = nothing
        #"--output"
        #    default=nothing
    end
    parsed_args = parse_args(args, s)

    defpars = Dict{String,Any}([
        "N"     => Dict("value" => 200, "type" => Int64),
        "k"     => Dict("value" => 4, "type" => Int64),
        "r"     => Dict("value" => 1.5,     "type" => Float64),
        "z" => Dict("value" => 1.3, "type" => Float64),
        "d" => Dict("value" => 0.9, "type" => Float64),
        "w" => Dict("value" => 0.01, "type" => Float64),
        "good_reputations" => Dict("value" => "payback", "type" => String),
        "num_trials" => Dict("value" => 10, "type" => Int64),
        "num_gens" => Dict("value" => 200, "type" => Int64),
        "gen_cutoff" => Dict("value" => 100, "type" => Int64),
        "output" => Dict("value" => "output/test.csv", "type" => String)
    ])
    pars = read_parameters(defpars, parsed_args["input"])

    # take the Cartesian product of all parameter combinations
    parsets = collect(Base.product(values(pars)...))
    nsets = length(parsets)

    # setup workers assuming directory is manually added to LOAD_PATH
    addprocs(min(parsed_args["ncpus"], round(Int64, Sys.CPU_THREADS / 2)))
    wpool = WorkerPool(workers())
    #extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)[1]
    extradir = filter((p)->match(r"/", p) !== nothing, LOAD_PATH)
    #@everywhere workers() push!(LOAD_PATH, $extradir)
    [@everywhere workers() push!(LOAD_PATH, $x) for x in extradir]
    @everywhere workers() eval(:(using Random))
    @everywhere workers() eval(:(using Statistics))
    @everywhere workers() eval(:(using NetworkLending))
    @everywhere workers() eval(:(using Dates))

    inputs  = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))
    results = RemoteChannel(()->Channel{Dict}(2 * nsets * maximum(pars["num_trials"])))

    @everywhere function run_worker(inputs, results)
        # save trial number and random seed
        seed = Dict(zip(["seed1", "seed2", "seed3", "seed4"], Random.GLOBAL_RNG.seed))

        while true
            pard = take!(inputs)
            pard = merge(pard, seed)

            N = pard["N"]
            k = pard["k"]
            r = pard["r"]
            z = pard["z"]
            d = pard["d"]
            w = pard["w"]
            if pard["good_reputations"] == "payback"
                good_reputations = [1,3]
            elseif pard["good_reputations"] == "coop"
                good_reputations = [2,3]
            elseif pard["good_reputations"] == "any"
                good_reputations = [1,2,3]
            elseif pard["good_reputations"] == "both"
                good_repuations = [3]
            end

            num_trials = pard["num_trials"]
            num_gens = pard["num_gens"]
            gen_cutoff = pard["gen_cutoff"]
            output = pard["output"]

            println("--- running ", pard["nrun"], " --- ")
            flush(stdout)

            gp = GameParams(r, z, d)
            game = NetworkGame(w, gp, good_reputations)
            network = Network(N, k)

            # initialize the NetworkPopulation object
            pop = NetworkPopulation(network, game)

            # initialize the frequency trajectories
            ft = FreqTraj(pop)

            for i in 1:num_gens
            	evolve_and_track!(pop, ft)
            end

            mean_freqs = mean(ft.freqs[end-gen_cutoff:end,:], dims=1)
            #mean_freqs = [sum(indv==x for indv in pop.lattice)*1.0/N^2 for x in types_to_compete]
            pard["mean_freqs"] = mean_freqs

            # return data to master process
            put!(results, pard)
        end
    end

    total_time_start = now()

    # load parameter sets into inputs channel
    nruns = 0
    for parset in parsets
        pard = Dict(zip(keys(pars), parset))
        #println(pard)
        println("--- queueing --- ")
        foreach(k->print(k, ": ", pard[k], ", "), sort(collect(keys(pard))))
        println()
        flush(stdout)
        for rep in 1:pard["num_trials"]
            nruns += 1
            rpard = copy(pard)
            rpard["rep"] = rep
            rpard["nrun"] = nruns
            put!(inputs, rpard)
        end
    end

    # start workers running on parameter sets in inputs
    for w in workers() # start tasks on the workers to process requests in parallel
        remote_do(run_worker, w, inputs, results)
    end

    # create output file name and data table
    output = pars["output"][1]
    println(output)
    file = occursin(r"\.csv$", output) ? output : output * ".csv"
    cols = push!(sort(collect(keys(pars))),
                 ["rep", "mean_freqs", "seed1", "seed2", "seed3", "seed4"]...)
    dat = DataFrame(Dict([(c, Any[]) for c in cols]))

    # grab results and output to CSV
    for sim in 1:nruns
        # get results from parallel jobs
        flush(stdout)
        resd = take!(results)
        nrun = pop!(resd, "nrun")

        # add to table (must convert dict keys to symbols) and save
        push!(dat, Dict([(Symbol(k), resd[k]) for k in keys(resd)]))
        CSV.write(file, dat)
    end
    total_time_stop = now()
    total_time = Dates.canonicalize(Dates.CompoundPeriod(round(total_time_stop - total_time_start, Dates.Second(1))))
    println("total time elapsed: $total_time")
end

#main(ARGS)

main(["--input", "submit/test_r_z_big.json"])
