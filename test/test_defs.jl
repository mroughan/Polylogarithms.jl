# common bits that the tests use
using Polylogarithms
using SpecialFunctions
using Test
using DataFrames, CSV
using Printf

import Base.MathConstants: π, pi, ℯ, e, γ, eulergamma, catalan, φ, golden

# functions for comparisons of results

# ≈   \approx
# x ≉ y is equivalent to !isapprox(x,y)
# isapprox(x, y; rtol::Real=atol>0 ? 0 : √eps, atol::Real=0, nans::Bool=false, norm::Function)

# but I want my own defs where I know what they are doing


# comon definitions for testing
machine_precision = eps(Float64)
accuracy_goal1 = 1e-12 # default low-precision goal
accuracy_goal2 = 1e-15 # default high-precision goal
relerr(z, x) = z == x ? 0.0 : abs(z - x) / abs(x)
relerrc(z, x) = max(relerr(real(z),real(x)), relerr(imag(z),imag(x)))
function near_equal(x, z, accuracy_goal)
    if typeof(x) <: AbstractFloat || typeof(z) <: AbstractFloat || typeof(x) <: Complex || typeof(z) <: Complex
        # for most cases we want small relative errors
        # but near zero, you run the risk of (1.0e-20 - 0.0)/0.0, or (0.0 - 1.0e-20)/1.0e-20
        # also, relative error in (a-b) might be bigger than the actual relative error in the terms
        if abs(z) < 1.0e-4 && abs(x - z) ≤ accuracy_goal 
            return true
        else
            return relerr(x,z) ≤ accuracy_goal
        end
    elseif typeof(x) <: Number && typeof(z) <: Number # deal with integers and rationals
        return x==z
    else
        # other things like tuples
        return x==z
    end
end

≈(a,b) = near_equal(a,b,accuracy_goal1)
≅(a,b) = near_equal(a,b,accuracy_goal2)
≒(a,b) = a === b
