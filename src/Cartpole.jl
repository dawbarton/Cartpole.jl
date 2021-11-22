module Cartpole

using OrdinaryDiffEq: ODEProblem, solve, Tsit5  # ODE integrator
using StaticArrays: SVector  # Fast statically-sized arrays
using MatrixEquations: arec  # Solve a Riccati equation
using Plots: plot  # Basic plotting
using ForwardDiff: jacobian  # Automatic differentiation
using LinearAlgebra: diagm, I  # Basic linear algebra

export cartpole, cartpole_defaults, cartpole_example

"""
    cartpole(u, p, t)

Implement the cart-pole equations of motion where `u = [x, x′, θ, θ′]` with `M` as the cart
mass, `m` as the pendulum mass, `L` as the length of the pendulum, and `g` as gravity.

```math
(M + m)x″ + mLsin(θ)θ′² - mLcos(θ)θ″ = f(u, p, t)
```

```math
mL²θ″ - mLcos(θ)x″ - mgLsin(θ) = 0
```

`f(u, p, t)` is an arbitrary control force.
"""
function cartpole(u, p, t)
    # Unpack the state
    (x, x′, θ, θ′) = u
    # Force to apply (a function)
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

function cartpole_integrate(u0, p, tend)
    # This assume Float64 as the underlying floating point type - could write this
    # generically which might be better for GPU implementations
    prob = ODEProblem(cartpole, SVector{4}(u0), (0.0, tend), p)
    return solve(prob, Tsit5())
end

"""
    cartpole_optimal_control(u, p, Q, R)

Calculate the optimal control gains (LQR) around a specific operating point.
"""
function cartpole_optimal_control(u, p, Q, R)
    A = jacobian(uu -> cartpole(uu, p, false), u)
    B = jacobian(pp -> cartpole(u, merge(p, (f = (args...) -> only(pp),)), false), [0.0])  # this is slightly cumbersome since I've made f as a function
    # Solve the Riccati equation
    P, _ = arec(A, B, R, Q)  # returns a tuple; the first element is the one of interest
    # Calculate optimal gains
    K = R \ (B' * P)
    return K
end

"""
    cartpole_example(; plotting)

Calculate the solution of the cart-pole example with optimal (LQR) control. Plotting is
optional, use `cartpole_example(plotting=true)` to display.
"""
function cartpole_example(; plotting = false)
    u0 = [0, 0, 0.1, 0]
    p0 = cartpole_defaults()
    tend = 10.0
    R = Matrix(I, 1, 1)
    Q = diagm([500.0, 250, 1, 1])
    K = cartpole_optimal_control(u0, p0, Q, R)  # calculate control gains
    p = merge(p0, (f = (u, p, t) -> only(-p.K * u), K = K))  # add a linear controller
    prob = ODEProblem(cartpole, SVector{4}(u0), (0.0, tend), p)
    sol = solve(prob, Tsit5())
    plotting && display(plot(sol))
    return sol
end

end  # module
