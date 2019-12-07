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
	)
	a = [0 r/2;
	r/2-1 r-1]
	return a
end

function b_matrix(
	a::Array{Float64, 2},
	k::Int64
	)
	b = zeros(Float64, size(a))
	for i in 1:size(b)[1]
		for j in 1:size(b)[2]
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
	function NOParams(r::Float64, k::Int64)
		a = a_matrix(r)
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
	[dx[i] = (x .* ((f + g) .- phi))[i] for i in 1:2]
	#println(dx)
end

k_vals = [3, 4, 5, 6]
w = 0.01 # selection strength
r_vals = collect(1:0.01:2)
coop_freq = zeros(Float64, length(k_vals), length(r_vals))
for (ki, k) in enumerate(k_vals)
	for (ri, r) in enumerate(r_vals)
		p = NOParams(r, k)
		x0 = [0.5, 0.5]
		tspan = (0.0,10000.0)
		prob = ODEProblem(dot_x!, x0, tspan, p)
		sol = solve(prob)
		coop_freq[ki, ri] = sol.u[end][2]
	end
end

fig = plt.figure()
[plt.plot(r_vals, coop_freq[ki, :], label="k = $k") for (ki, k) in enumerate(k_vals)]
[plt.vlines(2*k/(k+1), 0, 1, ls="--") for (ki, k) in enumerate(k_vals)]
plt.legend(loc=2)
plt.title("long-term frequency dynamics")
plt.xlabel(L"r")
plt.ylabel(L"\rho_c")
fig.tight_layout()
display(fig)
plt.savefig("figures/replicator_no_lending.pdf")
