using Plots

include("library/ou2d.jl")
include("defaults.jl")

# Default parameters
θ = ou2d_default_pars()
T = ou2d_default_tvec()

# Number of SDE realisations
n = 30

# Initial conditions to examine
ic = [:conditioned, :fixed, :perturbed]
θ̂  = [ou2d_sample_otherpars(θ;ic=icᵢ) for icᵢ in ic]

## Row 1
row1 = begin
    # Plot simulations
    plts_sims = [plot(xlabel="Time [d]") for _ = 1:3]
    for i = eachindex(ic)
        sols = [ou2d_solve_sde(θ,T;ic=ic[i]) for _ = 1:n]
        [plot!(plts_sims[i],sols[j],idxs=2,c=:red,α=0.2,label="") for j = 1:n]
        [plot!(plts_sims[i],sols[j],idxs=1,c=:blue,α=0.2,label="") for j = 1:n]
    end

    # Plot distributions
    function plt_dist_coord(d)
        if std(d) == 0
            xplt = [0.0,1.0]
            yplt = mean(d) * ones(2)
        else
            yplt = range(quantile.(d,[0.001,0.999])...,201)
            xplt = pdf.(d,yplt)
            xplt = xplt / maximum(xplt)
        end
        return xplt,yplt
    end
    plts_dists = [plot(xaxis=:flip,xlim=(0.0,1.2),xticks=[],xlabel="Density") for _ = 1:3]

    for i = eachindex(ic)
        # Obtain and plot initial condition
        d_x, d_y = ou2d_ic(θ;ic=ic[i]).v
        plot!(plts_dists[i],plt_dist_coord(d_y)...,c=:red,
            frange=std(d_y) == 0  ? Inf : 0.0,fα=0.2,label="")
        plot!(plts_dists[i],plt_dist_coord(d_x)...,c=:blue,
            frange=std(d_x) == 0  ? Inf : 0.0,fα=0.2,label="")
    end

    plot(vcat([[plts_dists[i],plts_sims[i]] for i = 1:3]...)...,link=:y,box=:on,
        layout=@layout([a{0.1w} b c{0.1w} d e{0.1w} f]),size=(1000,280))
end

## Row 2
row2 = begin
    # Plot for each initial condition
    plts = [plot(xlabel="Time [d]") for _ = 1:3]
    for i = eachindex(ic)
        # Solve moment equations for each initial condition
        sol1 = ou2d_solve_ode(θ,T;ic=ic[i])
        sol2 = ou2d_solve_ode(θ̂[i],T;ic=ic[i])
        # Plot
        plot!(plts[i],sol1,idxs=[1,3],c=[:blue :orange],lw=2.0,label=["Mean" "Variance"])
        plot!(plts[i],sol2,idxs=[1,3],c=:black,ls=:dash,lw=2.0,label="")
    end
    # Create row
    plot(plts...,link=:all,box=:on,layout=grid(1,3),size=(1000,280))
end

fig2 = plot(row1,row2,layout=grid(2,1),size=(1000,500))

savefig("$(@__DIR__)/fig2.svg")

plot!(fig2,xformatter=_->"",yformatter=_->"")

savefig("$(@__DIR__)/fig2_blank.svg")