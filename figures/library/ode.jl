using DifferentialEquations

include("shared.jl")

function ou2d_rhs!(dx,x,θ,t)
    a,b,c,d,e,f,p,r,s = θ
    m₁₀,m₀₁,m₂₀,m₀₂,m₁₁ = x
    dx[1] = a * e + b * f - a * m₁₀ - b * m₀₁
    dx[2] = c * e + d * f - c * m₁₀ - d * m₀₁
    dx[3] = 2((a * e + b * f) * m₁₀ - a * m₂₀ - b * m₁₁) + p^2
    dx[4] = 2((c * e + d * f) * m₀₁ - d * m₀₂ - c * m₁₁) + r^2 + s^2
    dx[5] = (c * e + d * f) * m₁₀ +  (a * e + b * f) * m₀₁ - c * m₂₀ - b * m₀₂ - (a + d) * m₁₁ + p * r
end

function ou2d_solve_ode(θ,T;ic=:fixed)
    X₀ = ou2d_ic(θ;ic,raw=true)
    solve(ODEProblem(ou2d_rhs!,X₀,extrema(T),θ),saveat=T)
end