using Plots
using Random

include("library/ou2d.jl")
include("defaults.jl")

Random.seed!(48)

# Parameters and time vector
θ = ou2d_default_pars()
T = range(0.0,10.0,2001)

# Stationary initial condition (i.e., sample Gaussian process)
ic = :stationary

# Simulate the SDE
sol = [ou2d_solve_sde(θ,T;ic) for _ = 1:3]

# Find zero crossing times
function zero_crossing(sol)
    t = sol.t
    u = [sol(tᵢ)[1] for tᵢ in t]
    if u[1] < 0
        if !any(u .> 0.0)
            return Inf
        end
        idx = findfirst(u .> 0.0)
    else
        if !any(u .< 0.0)
            return Inf
        end
        idx = findfirst(u .< 0.0)
    end
    t₁,u₁ = t[idx-1],u[idx-1]
    t₂,u₂ = t[idx],u[idx]
    t₁ + u₂ * (t₁ - t₂) / (u₁ - u₂)
end
tcross = zero_crossing.(sol)

# Create figure 1

plt1 = hline([0.0],c=:black,lw=2.0,ls=:dash,legend=:none)
plt2 = hline([0.0],c=:black,lw=2.0,ls=:dash,legend=:none)

for i = eachindex(sol)
    
    # Figure 1: untranslated time
    plot!(plt1,sol[i],idxs=1,c=[:blue,:red,:orange][i],lw=1.5)
    scatter!(plt1,[tcross[i]],[0.0],c=[:blue,:red,:orange][i],msw=0.0,m=:diamond,ms=5)

    # Figure 1: translated time (conditioned)
    plot!(plt2,sol[i].t .- tcross[i],hcat(sol[i].u...)[1,:],c=[:blue,:red,:orange][i],lw=1.5)
    scatter!(plt2,[0.0],[0.0],c=:black,msw=0.0,m=:diamond,ms=5)    

end

plot!(plt1,xlabel="t",ylabel="x(t)")
plot!(plt2,xlabel="t̂",ylabel="x(t)")

fig1 = plot(plt1,plt2,xlim=(0.0,3.0),ylim=(-0.5,1.5),widen=true,size=(600,250))
add_plot_labels!(fig1)
savefig(fig1,"$(@__DIR__)/fig1.svg")
fig1