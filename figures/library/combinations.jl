using LinearAlgebra
using Optim

function ou2d_si_quantities(θ;ic=:conditioned)
    a,b,c,d,e,f,p,r,s,x₀ = θ
    if ic == :stationary
        return ou2d_si_quantities(θ;ic=:conditioned)
    elseif ic == :conditioned
        [x₀,a + d,a*d - b*c,e,p^2,(d*p - b*r)^2 + b^2*s^2]
    elseif ic == :perturbed
        [x₀,a,b*c,d,e,p,(d*p - b*r)^2 + b^2*s^2]
    elseif ic == :fixed
        y₀ = θ[11]
        [x₀, a + d,a*d - b*c,e,p^2,d*p - b*r,b^2*s^2,d*x₀ - b*y₀ + b*f - d*e]
    elseif ic == :independent
        [a + d,a*d - b*c,e,f,p^2,r^2 + s^2,(d*p - b*r)^2 + b^2*s^2,(c*p - a*r)^2 + a^2*s^2]
    end
end

function ou2d_sample_otherpars(θ;ic=:conditioned)
    P = ou2d_si_quantities(θ;ic)
    θ̂ = rand(11)
    discrepancy(θ̂) = norm(P - ou2d_si_quantities(θ̂;ic))
    res = optimize(discrepancy,θ̂;g_tol=1e-20)
    ε,θ̂ = res.minimum, res.minimizer
    ε > 1e-5 && error("Did not converge! Try again :( ")
    display("Converged with ε = $ε")
    return θ̂
end
