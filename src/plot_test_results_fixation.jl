#!/usr/bin/env julia

## plt_test_results.jl
##
## Author: Taylor Kessinger <tkess@sas.upenn.edu>
## Plot results from NetworkLending simulations.
## Look at dependence of type frequencies on r and z.

using CSV, PyPlot, Statistics, NetworkLending, JSON

# load simulation output as a dataframe
runs = CSV.read("output/test_lending_fixation.csv")

# dicts to store fixation probabilities
type_freqs = Dict{Tuple{Int64, Float64, Float64, Float64},Array{Float64, 1}}()
fixation_times = Dict{Tuple{Int64, Float64, Float64, Float64},Float64}()
c_freqs = Dict{Tuple{Int64, Float64, Float64, Float64},Float64}()
p_freqs = Dict{Tuple{Int64,  Float64, Float64, Float64},Float64}()


N = sort(unique(runs[:N]))[1]

# get unique values from the runs dataframe
z_vals = sort(unique(runs[:z]))
d_vals = sort(unique(runs[:d]))
r_vals = sort(unique(runs[:r]))
k_vals = sort(unique(runs[:k]))

param_combs = collect(Base.product(k_vals, z_vals, d_vals, r_vals))

for (pi, param_comb) in enumerate(param_combs)
    k, z, d, r = param_comb
    #println("$param_comb")
    type_freqs[param_comb] = zeros(Float64, 4)
    c_freqs[param_comb] = 0.0
    p_freqs[param_comb] = 0.0
    fixation_times[param_comb] = 0.0
    tmp_runs = runs[(runs[:k] .== k) .& (runs[:z] .== z) .& (runs[:d] .== d) .& (runs[:r] .== r), :]
    for (ri, run) in enumerate(eachrow(tmp_runs[:,:]))
        fixation_times[param_comb] += 1.0*run[:generations]/size(tmp_runs, 1)
        #println("$(run[:mean_freqs])")
        freqs = parse.(Float64,String.(split(run[:mean_freqs], r";|,| |\[|\]")[2:2:end-1]))
        #freqs = JSON.parse(run[:mean_freqs])
        #println("$freqs")
        type_freqs[param_comb] += freqs/size(tmp_runs, 1)
        c_freqs[param_comb] += sum([freqs[3], freqs[4]])/size(tmp_runs, 1)
        p_freqs[param_comb] += sum([freqs[2], freqs[4]])/size(tmp_runs, 1)
    end
end

zmin, zmax, rmin, rmax = z_vals[1], z_vals[end], r_vals[1], r_vals[end]
scaling = (zmax-zmin)/(rmax-rmin)

for (ki, k) in enumerate(k_vals)
    fig, axs = plt.subplots(4, length(d_vals), figsize=(12,12),
        sharex="col", sharey="row")
    for (di, d) in enumerate(d_vals)
        freqs_im = zeros(length(r_vals), length(z_vals), 4)
        for (zi, z) in enumerate(z_vals)
            for (ri, r) in enumerate(r_vals)
                freqs_im[ri, zi, :] += type_freqs[k, z, d, r]
            end
        end
        for i in 1:4
            ax = axs[i,di]
            im = ax.imshow(freqs_im[:,:,i], origin = "lower",
                aspect=scaling,
                vmin=0, vmax=1,
                extent = [zmin, zmax, rmin, rmax])
            if di == 1
                ax.set_ylabel("r")
            end
            if i == 1
                ax.set_title("d = $d")
            end
        end
        #fig.colorbar(im)
    #    cbar = ax.cax.colorbar(im)
    #    cbar = grid.cbar_axes[0].colorbar(im)
    end
    fig.suptitle("payback, k = $k")
    #plt.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #plt.colorbar(im, cax=cbar_ax)
    #fig.colorbar(im)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    #plt.subplots_adjust(top=0.85)
    display(fig)


    fig, axs = plt.subplots(2, length(d_vals), figsize=(9, 6),
        sharex="col", sharey="row")
    for (di, d) in enumerate(d_vals)
        c_freq_im = zeros(length(r_vals), length(z_vals))
        p_freq_im = zeros(length(r_vals), length(z_vals))
        for (zi, z) in enumerate(z_vals)
            for (ri, r) in enumerate(r_vals)
                c_freq_im[ri, zi] = c_freqs[k, z, d, r]
                p_freq_im[ri, zi] = p_freqs[k, z, d, r]
            end
        end
        ax = axs[1,di]
        im = ax.imshow(p_freq_im, origin = "lower",
            aspect=scaling,
            vmin=0, vmax=1,
            extent = [zmin, zmax, rmin, rmax])
        if di == 1
            ax.set_ylabel("r")
        end
        ax.set_title("d = $d")

        ax = axs[2,di]
        im = ax.imshow(c_freq_im, origin = "lower",
            aspect=scaling,
            vmin=0, vmax=1,
            extent = [zmin, zmax, rmin, rmax])
        ax.set_xlabel("z")
        if di == 1
            ax.set_ylabel("r")
        end
        #fig.colorbar(im)
    #    cbar = ax.cax.colorbar(im)
    #    cbar = grid.cbar_axes[0].colorbar(im)
    end
    fig.suptitle("payback, k = $k")
    #plt.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #plt.colorbar(im, cax=cbar_ax)
    #fig.colorbar(im)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    #plt.subplots_adjust(top=0.85)
    display(fig)



    fig, axs = plt.subplots(1, length(d_vals), figsize=(9, 6),
        sharex="col", sharey="row")
    for (di, d) in enumerate(d_vals)
        generations = zeros(length(r_vals), length(z_vals))
        for (zi, z) in enumerate(z_vals)
            for (ri, r) in enumerate(r_vals)
                generations[ri, zi] = fixation_times[k, z, d, r]
            end
        end
        ax = axs[di]
        im = ax.imshow(generations, origin = "lower",
            aspect=scaling,
            extent = [zmin, zmax, rmin, rmax])
        if di == 1
            ax.set_ylabel("r")
        end
        ax.set_xlabel("z")
        ax.set_title("d = $d")
        fig.colorbar(im, ax=ax)
    end
    fig.suptitle("payback, k = $k")
    #plt.subplots_adjust(right=0.8)
    #cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    #plt.colorbar(im, cax=cbar_ax)
    #fig.colorbar(im)
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    #plt.subplots_adjust(top=0.85)
    display(fig)
end
