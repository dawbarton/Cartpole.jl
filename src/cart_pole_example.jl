module cart_pole_example

using OrdinaryDiffEq  # ODE integrator
using StaticArrays  # Fast statically-sized arrays
using MatrixEquations  # Solve a Riccati equation
using Plots  # Basic plotting
using ForwardDiff  # Automatic differentiation

export cartpole, cartpole_defaults, cartpole_example

"""
    cartpole(u, p, t)

Implement the cart-pole equations of motion where `u = [x, x′, θ, θ′]` with `M` as the cart
mass, `m` as the pendulum mass, `L` as the length of the pendulum, and `g` as gravity.

```math
(M + m)x″ + mLsin(θ)θ′² - mLcos(θ)θ″ = f(u, p, t)
mL²θ″ - mLcos(θ)x″ - mgLsin(θ) = 0
```

`f(u, p, t)` is an arbitrary control force.
"""
function cartpole(u, p, t)
    # Unpack the state
    (x, x′, θ, θ′) = u
    # Force to apply
    force = p.f(u, p, t)
    # Right-hand side
    return SVector(
        x′,
        (force + p.m * p.L * sin(θ) * θ′^2 + p.m * p.g * cos(θ) * sin(θ)) /
        (p.M + p.m - p.m * cos(θ)^2),
        θ′,
        (
            force * cos(θ) -
            (p.M + p.m) * p.g * sin(θ) +
            p.m * p.L * cos(θ) * sin(θ) * x′^2
        ) / (p.m * p.L * cos(θ)^2 - (p.M + p.m) * p.L),
    )
end

function cartpole_defaults()
    return (M = 5, m = 1, L = 2, g = 10, f = (u, p, t) -> false)
end

function cartpole_example(u0, tend, p = cartpole_defaults())
    prob = ODEProblem(cartpole, SVector{4}(u0), (0.0, tend), p)
    return solve(prob, Tsit5())
end

end  # module
