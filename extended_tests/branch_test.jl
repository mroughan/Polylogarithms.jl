using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot
using Printf
include("utilities.jl")

series_1 = Polylogarithms.polylog_series_1
series_2 = Polylogarithms.polylog_series_2
series_3 = Polylogarithms.polylog_series_3
polylog_reciprocal = Polylogarithms.polylog_reciprocal 

# test behaviour along the branch
s = -1.5
z = 3.0 .+ collect( -0.01: 0.001: 0.01 )*im

for i=1:length(z)
    println(" ", imag(z[i]), ", ", imag( series_2(s, z[i])[1]), ", ", imag( polylog_reciprocal(s, z[i])[1]) )
end




