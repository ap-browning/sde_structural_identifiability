using Plots

include("library/ou2d.jl")
include("defaults.jl")

# Parameters and time vector
θ = ou2d_default_pars()
θ̂ = ou2d_sample_otherpars(θ;ic=:independent)

ρ = ou2d_acov(θ)
ρ̂ = ou2d_acov(θ̂)

fig3 = plot()

plot!(fig3, t -> ρ(t)[1,1],xlim=(0.0,22.0),lw=2.0,c=:blue)
plot!(fig3, t -> ρ̂(t)[1,1],xlim=(0.0,22.0),lw=2.0,c=:black,ls=:dash)

plot!(fig3, t -> ρ(t)[2,2],xlim=(0.0,22.0),lw=2.0,c=:red)
plot!(fig3, t -> ρ̂(t)[2,2],xlim=(0.0,22.0),lw=2.0,c=:black,ls=:dash)

plot!(fig3, t -> ρ(t)[1,2],xlim=(0.0,22.0),lw=2.0,c=:orange)
plot!(fig3, t -> ρ̂(t)[1,2],xlim=(0.0,22.0),lw=2.0,c=:orange,ls=:dash)

plot!(fig3,xlim=(0.0,20.0),widen=true)
plot!(fig3,xlabel="Time [h]")
plot!(fig3,size=(300,250))

savefig(fig3,"$(@__DIR__)/fig3.svg")
fig3