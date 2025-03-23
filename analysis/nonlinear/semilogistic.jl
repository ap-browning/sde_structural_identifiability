using StructuralIdentifiability

# Equivalent ODE
ode = @ODEmodel(
    x1'(t) = a * x1(t) * (1 - b * x1(t)) + c * x2(t),
    x2'(t) = d * x1(t) * (1 - e * x1(t)) + f * x2(t),
    y1(t) = x1(t)
)
assess_identifiability(ode)
res = find_identifiable_functions(ode, with_states = true)