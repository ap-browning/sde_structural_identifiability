using StructuralIdentifiability

# Equivalent ODE
ode = @ODEmodel(
    x1'(t) = a * x1(t) + b * x2(t),
    x2'(t) = c * x2(t) + d * x1(t) * x2(t),
    y1(t) = x1(t)
)
assess_identifiability(ode)
res = find_identifiable_functions(ode, with_states = true)