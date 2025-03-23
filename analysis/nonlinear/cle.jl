using StructuralIdentifiability

# Equivalent ODE
ode = @ODEmodel(
    x1'(t) = -2α * x1(t)^2 + 2β * x2(t) + ε - ζ * x1(t),
    x2'(t) = α * x1(t)^2 + γ - (β + δ) * x2(t),
    y1(t) = x1(t)
)
assess_identifiability(ode)
res = find_identifiable_functions(ode, with_states = true)