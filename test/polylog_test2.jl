# test using data
using Polylogarithms
using SpecialFunctions
using Test
using DataFrames, CSV
import Base.MathConstants: π, pi, ℯ, e, γ, eulergamma, catalan, φ, golden
include("test_defs.jl")
include("../bench/utilities.jl")


# desired accuracy is 1.0e-12, but we get a few points above this so
≈(a,b) = near_equal(a,b,1.0e-10)
