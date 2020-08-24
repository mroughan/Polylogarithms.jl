"""
    Polylogarithms

Module containing functions to calculate the polylogarithm and associated functions

"""
module Polylogarithms
import SpecialFunctions

# not using MPFR for the moment
# using Base.MPFR: ROUNDING_MODE, big_ln2

export polylog, bernoulli, harmonic, stieltjes, dirichlet_beta
export Diagnostics

# @compat ComplexOrReal{T} = Union{T,Complex{T}}
# s::ComplexOrReal{Float64}
ComplexOrReal{T} = Union{T,Complex{T}}

# Constants
include("constants.jl")

# Series
include("stieltjes.jl")
include("bernoulli_n.jl")
include("harmonic.jl")

# Functions
include("beta.jl")
include("bernoulli_poly.jl")
include("polylog.jl")

end
