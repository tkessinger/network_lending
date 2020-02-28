#!/usr/bin/env julia

## analytical_payoffs.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Analytically solve the replicator equation for our network lending model.

using Revise
using NetworkLending
using PyPlot
using DifferentialEquations

function a_matrix(
	r::Float64,
	z::Float64,
	d::Float64
	)
	a = zeros(Float64, 4, 4)
	gp = GameParams(r, z, d)
	for i in 0:3
		for j in 0:3
			a[i+1,j+1] = pair_PGG_payoffs([i, j], [i%2, j%2], gp)[1]
		end
	end
	return a
end

function a_matrix_long_version(
	r::Float64,
	z::Float64,
	d::Float64
	)
	p = [0 0 1-d 1;
		0 0 1-d 1;
		1-d 1-d 2-2d 2-d;
		1 1 2-d 2]
	l = [0 0 0 0;
		1 1 1 1;
		0 0 0 0;
		1 1 1 1]
	R = [0 0 0 0;
		0 0 0 0;
		1-d 1-d 1-d 1-d;
		1 1 1 1]
	a = (r/2*p - R - z*l)
	return a
end

function a_matrix_new_version(
	r::Float64,
	z::Float64,
	d::Float64
	)
	a = [1-d 1-d (1-d)*(1+r/2) 1-d+r/2;
	-z -z r/2*(1-d)-z r/2-z;
	r*(1-d)/2 r*(1-d)/2 r*(1-d) (r/2)*(2-d);
	r/2-1-z r/2-1-z r/2*(2-d)-1-z r-1-z]
	# p = [0 0 1-d 1;
	# 	0 0 1-d 1;
	# 	1-d 1-d 2-2d 2-d;
	# 	1 1 2-d 2]
	# l = [0 0 0 0;
	# 	1 1 1 1;
	# 	0 0 0 0;
	# 	1 1 1 1]
	# R = [0 0 0 0;
	# 	0 0 0 0;
	# 	1-d 1-d 1-d 1-d;
	# 	1 1 1 1]
	# a = (r/2*p - R - z*l)
	return a
end

function b_matrix(
	a::Array{Float64, 2},
	k::Int64
	)
	b = zeros(Float64, 4, 4)
	for i in 1:4
		for j in 1:4
			b[i, j] = 1.0/((k+1)*(k-2))*((k+1)*(a[i,i] - a[j,j]) + a[i,j] - a[j,i])
		end
	end
	return b
end

struct NOParams
	# nuisance structure for holding the payoff matrix
	# and Nowak-Ohtsuki transform
	a::Array{Float64, 2}
	b::Array{Float64, 2}

	# constructors
	function NOParams(gp::GameParams, k::Int64)
		a = a_matrix(gp.r, gp.z, gp.d)
		b = b_matrix(a, k)
		return new(a, b)
	end
end

function dot_x!(
	dx::Array{Float64, 1},
	x::Array{Float64, 1},
	p::NOParams,
	t::Float64
	)
	a, b = p.a, p.b
	f = a*x
	g = b*x
	phi = x' * f
	[dx[i] = (x .* ((f + g) .- phi))[i] for i in 1:4]
	#println(dx)
end

k = 4 # network degree
w = 0.01 # selection strength
r_vals = collect(1:0.05:2)
z_vals = collect(0:0.05:2)
d_vals = [0.0, 0.5, 0.9, 0.99]

zmin, zmax = minimum(z_vals), maximum(z_vals)
rmin, rmax = minimum(r_vals), maximum(r_vals)
scaling = (zmax-zmin)/(rmax-rmin)

initial_types = [2, 3, 0]

strats = join(initial_types)

fig, axs = plt.subplots(4, length(d_vals), figsize=(12,12),
	sharex="col", sharey="row")

save_fig = true
make_fig = true

	if make_fig
	for (di, d) in enumerate(d_vals)
		results = zeros(Float64, length(r_vals), length(z_vals), 4)
		for (ri, r) in enumerate(r_vals)
			for (zi, z) in enumerate(z_vals)
				println("$r, $z")
				gp = GameParams(r, z, d)
				p = NOParams(gp, k)
				x0 = zeros(Float64, 4)
				[x0[i+1] = 1.0/length(initial_types) for i in initial_types]
				tspan = (0.0,10000.0)
				prob = ODEProblem(dot_x!, x0, tspan, p)
				sol = solve(prob)
				results[ri, zi, :] = sol.u[end]
			end
		end
		for i in 1:4
			ax = axs[i, di]
			ax.imshow(results[:,:, i], origin="lower", vmin=0, vmax=1,
			 	extent = [zmin, zmax, rmin, rmax], aspect = scaling)
			if di == 1
				ax.set_ylabel("r")
			end
			if i == 1
				ax.set_title("d = $d")
			end
			if i == 4
				ax.set_xlabel("z")
			end
		end
	end
	fig.suptitle("numerical results, strats = $strats, k = $k")
	#plt.subplots_adjust(right=0.8)
	#cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
	#plt.colorbar(im, cax=cbar_ax)
	#fig.colorbar(im)
	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
	display(fig)
	if save_fig
		plt.savefig("figures/new_numerical_results_k_$(k)_strats_$strats.pdf")
	end
end
# for (di, d) in enumerate(d_vals)
# 	payback_freq = zeros(Float64, length(r_vals), length(z_vals))
# 	coop_freq = zeros(Float64, length(r_vals), length(z_vals))
# 	results = zeros(Float64, length(r_vals), length(z_vals), 4)
# 	for (ri, r) in enumerate(r_vals)
# 		for (zi, z) in enumerate(z_vals)
# 			println("$r, $z")
# 			gp = GameParams(r, z, d)
# 			p = NOParams(gp, k)
# 			x0 = zeros(Float64, 4)
# 			[x0[i+1] = 1.0/length(initial_types) for i in initial_types]
# 			tspan = (0.0,10000.0)
# 			prob = ODEProblem(dot_x!, x0, tspan, p)
# 			sol = solve(prob)
# 			results[ri, zi, :] = sol.u[end]
# 			payback_freq[ri, zi] = sol.u[end][2] + sol.u[end][4]
# 			coop_freq[ri, zi] = sol.u[end][3] + sol.u[end][4]
# 		end
# 	end
#
# 	zmin, zmax, rmin, rmax = z_vals[1], z_vals[end], r_vals[1], r_vals[end]
# 	scaling = (zmax-zmin)/(rmax-rmin)
#
# 	savefig = true

	# fig = plt.figure()
	# v
	# plt.xlabel(L"z")
	# plt.ylabel(L"r")
	# plt.title("coop frequency, k = $k, d = $d, initial types = $strat_string")
	# plt.colorbar()
	# plt.tight_layout()
	# display(fig)
	# if savefig
	# 	plt.savefig("figures/new_coop_frequency_k_$(k)_d_$(d)_strats_$strat_string.pdf")
	# end
	#
	# fig = plt.figure()
	# plt.imshow(payback_freq, origin="lower", vmin=0, vmax=1,
	# 	extent = [zmin, zmax, rmin, rmax], aspect = scaling)
	# plt.xlabel(L"z")
	# plt.ylabel(L"r")
	# plt.title("payback frequency, k = $k, d = $d, initial types = $strat_string")
	# plt.colorbar()
	# plt.tight_layout()
	# display(fig)
	# if savefig
	# 	plt.savefig("figures/new_payback_frequency_k_$(k)_d_$(d)_strats_$strat_string.pdf")
	# end

# 	fig, axs = plt.subplots(2,2, sharey="col", sharex="row",
# 		figsize = (8,8))
#
# 	for i in 1:4
# 		ax = axs[i]
# 		im = ax.imshow(results[:,:,i], origin="lower", vmin=0, vmax=1,
# 			extent = [zmin, zmax, rmin, rmax], aspect=scaling)
# 		if i ∈ [2, 4]
# 			ax.set_xlabel(L"z")
# 		end
# 		if i ∈ [1, 2]
# 			ax.set_ylabel(L"r")
# 		end
# 		ax.set_title("frequency of type $(i-1)")
# 	end
# 	fig.suptitle("type frequencies, k = $k, d = $d, initial types = $strat_string")
# 	fig.tight_layout(rect=[0, 0.03, 1, 0.96])
# 	display(fig)
# 	if savefig
# 		plt.savefig("figures/analytic_type_freqs_k_$(k)_strats_$strat_string.pdf")
# 	end
#
# end
#
# fig = plt.figure(figsize = (5,5))
# for (di, d) in enumerate(d_vals)
# 	predicted_r_vals = 2*k .* (1.0 .+ z_vals) ./ (k+1.0-d)
# 	plt.plot(z_vals, predicted_r_vals, label="d = $d")
# end
# plt.xlabel(L"z")
# plt.ylabel(L"r")
# plt.legend(loc=2)
# plt.title("theoretical prediction, k = $k")
# plt.tight_layout()
# display(fig)
# if savefig
# 	plt.savefig("figures/r_prediction_k_$(k)_strats_$strat_string.pdf")
# end
