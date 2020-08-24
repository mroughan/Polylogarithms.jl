# test using data
using Polylogarithms
using SpecialFunctions
using Test
using DataFrames, CSV
import Base.MathConstants: π, pi, ℯ, e, γ, eulergamma, catalan, φ, golden
include("test_defs.jl")
using Printf

# desired accuracy is 1.0e-12, but we get a few points above this so
≈(a,b) = near_equal(a,b,1.0e-10)

L = Symbol("Li_s(z)")
@testset "Polylogarithm polylog on data" begin
    for C=1:3
        filename = @sprintf("../data/polylog_test_data_rand_%d.csv", C)
        data1 = CSV.read(filename, DataFrame; delim=",", type=String)

        # has trouble reading in numbers like "2." so read all into strings, and parse
        data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
        data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
        data1[!, L] = parse.(Complex{Float64}, data1[!, L] )
        
        m = size(data1,1)
        Li = data1[!, L]
        s  = data1[!,:s]
        z  = data1[!,:z]

        for i=1:m
            @test polylog( s[i], z[i] ) ≈ Li[i]
        end

    end
end
