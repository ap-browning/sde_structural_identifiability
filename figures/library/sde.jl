using DifferentialEquations

include("shared.jl")

function ou2d_f!(dx,x,Θ,t)
    A,B,S = Θ
    dx .= -A * (x - B)
end
function ou2d_g!(dx,x,Θ,t)
    A,B,S = Θ
    dx .= S
end

function ou2d_solve_sde(θ,T;ic=:fixed)
    X₀ = ou2d_ic(θ;ic)
    Θ = ou2d_pars_to_mat(θ)
    prob = SDEProblem(ou2d_f!,ou2d_g!,rand(X₀),extrema(T),Θ,noise_rate_prototype=zeros(2,2))
    solve(prob,saveat=T,tstops=T)
end

θ = ou2d_default_pars()
T = ou2d_default_tvec()

sol = ou2d_solve_sde(θ,T)