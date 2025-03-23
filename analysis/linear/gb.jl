using StructuralIdentifiability

## O(1)
ode = @ODEmodel(
    x1'(t) = a * e + b * f - a * x1(t) - b * x2(t),
    x2'(t) = c * e + d * f - c * x1(t) - d * x2(t),
    y1(t) = x1(t)
)
assess_identifiability(ode)
res = find_identifiable_functions(ode, with_states = true)

## O(2)
ode = @ODEmodel(
    x1'(t) = a * e + b * f - a * x1(t) - b * x2(t),
    x2'(t) = c * e + d * f - c * x1(t) - d * x2(t),
    x3'(t) = 2((a*e + b*f) * x1(t) - a * x3(t) - b * x5(t)) + p^2 * x3(t),
    x4'(t) = 2((c*e + d*f) * x2(t) - d * x4(t) - c * x5(t)) + (r^2 + s^2) * x4(t),
    x5'(t) = (a*e + b*f)*x2(t) - b * x4(t) + (c*e + d*f)*x1(t) - (a + d - p*r) * x5(t) - c*x3(t),
    y1(t) = x1(t),
    y2(t) = x3(t)
)
assess_identifiability(ode)
res = find_identifiable_functions(ode, known_ic = [e], with_states = true)

## O(3)
ode = @ODEmodel(
    x1'(t) = a * e + b * f - a * x1(t) - b * x2(t),
    x2'(t) = c * e + d * f - c * x1(t) - d * x2(t),
    x3'(t) = 2((a*e + b*f) * x1(t) - a * x3(t) - b * x5(t)) + p^2 * x3(t),
    x4'(t) = 2((c*e + d*f) * x2(t) - d * x4(t) - c * x5(t)) + (r^2 + s^2) * x4(t),
    x5'(t) = (a*e + b*f)*x2(t) - b * x4(t) + (c*e + d*f)*x1(t) - (a + d - p*r) * x5(t) - c*x3(t),
    x6'(t) = 3*((a*e + b*f) * x3(t) + (p^2 - a) * x6(t) - b * x8(t)),
    x7'(t) = 3*((c*e + d*f) * x4(t) + (r^2 + s^2 - d)*x7(t) - c * x9(t)),
    x8'(t) = (c*e + d*f)*x3(t) + 2*(a*e + b*f)*x5(t) - c*x6(t) + (-2a - d + p^2 + 2p*r)*x8(t) - 2b * x9(t),
    x9'(t) = (a*e + b*f)*x4(t) + 2*(c*e + d*f)*x5(t) - b*x7(t) -  2c*x8(t) + (-a - 2d + 2p*r + r^2 + s^2)*x9(t),
    y1(t) = x1(t),
    y2(t) = x3(t),
    y3(t) = x6(t)
)
assess_identifiability(ode)
res = find_identifiable_functions(ode, with_states = true)