#!/usr/bin/env julia

## test_network_lending.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Test NetworkLending and plot frequency trajectories.

using Revise
using NetworkLending
using PyPlot

N = 1000 # population size
k = 4 # network degree
w = 0.01 # selection strength

good_reputations = [1, 3]

r_vals = [1.2, 1.4, 1.6, 1.8]
z_vals = [0.4, 0.8, 1.2, 1.6]
d_vals = [0.0, 0.5, 0.9, 0.99]

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

colors = ["k", "y", "b", "g"]

for (di, d) in enumerate(d_vals)
	fig, axs = plt.subplots(length(r_vals), length(z_vals),
		sharey = "row", figsize=(12,12))
	for (ri, r) in enumerate(r_vals)
		for (zi, z) in enumerate(z_vals)
			# initialize the game and network basics
			gp = GameParams(r, z, d)
			game = NetworkGame(w, gp, good_reputations)
			network = Network(N, k)

			# initialize the NetworkPopulation object
			pop = NetworkPopulation(network, game)

			# initialize the frequency trajectories
			ft = FreqTraj(pop)

			# for i in 1:10*N
			# 	evolve_and_track!(pop, ft)
			# end
			evolve_and_track_until_fixation!(pop, ft)
			ax = axs[ri, zi]
			[ax.plot(ft.freqs[:,x], label=("$(strategy_label(x-1))"), c = colors[x]) for x in 1:4]
			ax.set_title("r = $r, z = $z")
			if ri == length(r_vals)
				ax.set_xlabel("time")
			end
			if zi == length(z_vals)
				ax.set_ylabel(L"\rho")
			end
		end
	end
	plt.legend(loc=2)
	fig.suptitle("d = $d, k = $k")
	#plt.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	#plt.colorbar(im, cax=cbar_ax)
	#fig.colorbar(im)
	fig.tight_layout(rect=[0, 0, 1, 0.96])
	display(fig)
	plt.savefig("figures/traj_grid_d_$(d)_k_$(k).pdf")
end

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
