## 2D Ornstein-Uhlenbeck model
#
# dx = -A (x - B) dt + S dW
#
# for A = [a b; c d], B = [e; f] and S = [p 0; r s]

using Distributions
using LinearAlgebra
using Kronecker
using Random

## PARAMETER VECTOR TO MATRICES
function ou2d_pars_to_mat(θ)
    a,b,c,d,e,f,p,r,s = θ
    A = [a b; c d]
    B = [e; f]
    S = [p 0; r s]
    return A,B,S
end

## RANDOMLY CHOSEN PARAMETER SET
ou2d_default_pars() = rand(MersenneTwister(5),11)
ou2d_default_tvec() = range(0.0,2.0,201)

## INITIAL CONDITIONS
function ou2d_ic_stationary(A,B,S)
    if any(eigvals(A) .≤ 0.0)
        error("The system has no stationary distribution!")
    end
    Σ = reshape((A ⊕ A) \ (S * S')[:],2,2)
    return MvNormal(B[:],(Σ + Σ') / 2)
end
ou2d_ic_stationary(θ) = ou2d_ic_stationary(ou2d_pars_to_mat(θ)...)

function ou2d_ic_conditioned(dist::MvNormal,x₀)
    e,f = mean(dist)
    Σ = cov(dist)
    μ = f + Σ[1,2] / Σ[1,1] * (x₀ - e)
    σ = sqrt(Σ[2,2] - Σ[1,2]^2 / Σ[1,1])
    return Product([Normal(x₀,0.0),Normal(μ,σ)])
end
ou2d_ic_conditioned(A,B,S,x₀) = ou2d_ic_conditioned(ou2d_ic_stationary(A,B,S)...,x₀)
ou2d_ic_conditioned(θ::Vector,x₀) = ou2d_ic_conditioned(ou2d_ic_stationary(θ),x₀)
ou2d_ic_conditioned(θ::Vector) = ou2d_ic_conditioned(θ[1:9],θ[10])

function ou2d_ic_perturbed(d::MvNormal,x₀)
    d_x = Normal(x₀,0.0)
    d_y = Normal(mean(d)[2],sqrt(cov(d)[2,2]))
    return Product([d_x,d_y])
end
ou2d_ic_perturbed(A,B,S,x₀) = ou2d_ic_perturbed(ou2d_ic_stationary(A,B,S),x₀)
ou2d_ic_perturbed(θ::Vector,x₀) = ou2d_ic_perturbed(ou2d_ic_stationary(θ),x₀)
ou2d_ic_perturbed(θ::Vector) = ou2d_ic_perturbed(θ[1:9],θ[10])

function ou2d_ic_fixed(x₀,y₀)
    return Product([Normal(x₀,0.0),Normal(y₀,0.0)])
end
ou2d_ic_fixed(θ::Vector) = ou2d_ic_fixed(θ[10:11]...)

function raw_moments(d::Distribution)
    m₁₀,m₀₁ = mean(d)
    Σ = cov(d)
    m₂₀ = Σ[1,1] + m₁₀^2
    m₀₂ = Σ[2,2] + m₀₁^2
    m₁₁ = Σ[1,2] + m₁₀ * m₀₁
    return [m₁₀,m₀₁,m₂₀,m₀₂,m₁₁]
end

function ou2d_ic(θ::Vector;ic=:stationary,raw=false)
    if ic == :stationary
        d = ou2d_ic_stationary(θ)
    elseif ic == :conditioned
        d = ou2d_ic_conditioned(θ)
    elseif ic == :perturbed
        d = ou2d_ic_perturbed(θ)
    elseif ic == :fixed
        d = ou2d_ic_fixed(θ)
    end
    if raw
        return raw_moments(d)
    else
        return d
    end
end


function ou2d_acov(A,B,S)
    Σ∞ = reshape(inv(A ⊕ A) * (S * S')[:],2,2)
    t -> Σ∞ * exp(-A' * t)
end
ou2d_acov(θ) = ou2d_acov(ou2d_pars_to_mat(θ)...)