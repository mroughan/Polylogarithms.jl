"""
    Polylogarithms

Module containing functions to calculate the polylogarithm and associated functions

"""
module Polylogarithms
import SpecialFunctions
# not using MPFR for the moment
# using Base.MPFR: ROUNDING_MODE, big_ln2

export polylog, bernoulli, euler, harmonic, stieltjes, dirichlet_beta
export rogers, spence, dilog, trilog
export Diagnostics
export parse
export Memoization, CacheZeta, clearcache

# @compat ComplexOrReal{T} = Union{T,Complex{T}}
# s::ComplexOrReal{Float64}
ComplexOrReal{T} = Union{T,Complex{T}}

const default_accuracy = 1.0e-12
const default_max_iterations = 1000
const near_int_threshold = 1.0e-6
const series_transition_threshold = 0.25

# Some extra parsing routines for reading Mathematics output, but also, complex numbers
include("utilities.jl")

# Constants
include("constants.jl")

# Series
include("stieltjes.jl")
include("bernoulli_n.jl")
include("harmonic.jl")
# include("gamma_derivatives.jl") # this just generates a table that we include into the code

# Functions
include("beta.jl")
include("bernoulli_poly.jl")
include("polylog.jl")
include("additional_functions.jl")

end
