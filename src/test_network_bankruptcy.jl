#!/usr/bin/env julia

## test_network_bankruptcy.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test NetworkLending and plot frequency trajectories.

using Revise
using NetworkBankruptcy
using PyPlot

N = 10 # population size
k = 4 # network degree
r = 1.6 # synergy factor
d = 0.9 # decrement for bad actors
w = 0.05 # selection strength

# initialize the game and network basics
gp = GameParams(r, d)
game = NetworkGame(w, gp)
network = Network(N, k)

# initialize the NetworkPopulation object
pop = NetworkPopulation(network, game)
#pop.reputations = zeros(Int64, pop.network.N)

# initialize the frequency trajectories
ft = FreqTraj(pop)

# for i in 1:10*N
# 	evolve_and_track!(pop, ft)
# end

evolve_and_track_until_fixation!(pop, ft)

function strategy_label(
	strategy::Int64
	)
	return (strategy == 0 ? "defect" : "cooperate")
end

colors = []

fig = plt.figure()
[plt.plot(ft.freqs[:,x], label=("$(strategy_label(x-1))")) for x in 1:2]
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
