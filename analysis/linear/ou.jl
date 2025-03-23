using StructuralIdentifiability

## Observe only x(t)
ode = @ODEmodel(
    x1'(t) = a * e + b - a * x1(t) - f * x2(t),
    x2'(t) = c * e + d - c * x1(t) - d * x2(t),
    x3'(t) = p^2 + 2a * e * x1(t) + 2x1(t) - 2a * x3(t) - 2x5(t),
    x4'(t) = r^2 + s^2 + 2c * e * x2(t) + 2d * x2(t) - 2d * x4(t) - 2c * x5(t),
    x5'(t) = p * r + c * e * x1(t) + d * x1(t) + a * e * x2(t) + x2(t) - c * x3(t) - x4(t) - a * x5(t) - d * x5(t),
    y1(t) = x1(t),
    y2(t) = x3(t)
)
assess_identifiability(ode)
find_identifiable_functions(ode, with_states = true, known_ic = [p], simplify = :weak)