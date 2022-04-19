using DataFrames, CSV
using PyPlot
using Printf
using SpecialFunctions
include(joinpath(@__DIR__, "..", "src", "utilities.jl"))

filename = joinpath(@__DIR__, "..", "data", "zeta_test_data.csv")
data1 = CSV.read(filename; delim=",", type=String)

data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
data1[!,:zeta] = parse.(Complex{Float64}, data1[!,:zeta] )

m = size(data1,1)
s = data1[!,:s]
z = complex(data1[!,:z])

S1 = zeros(Complex{Float64}, m)
for i = 1:m
    S1[i] = SpecialFunctions.zeta( 1-s[i], 0.5 + log(z[i])/(2*pi*im) )
end

error1 = S1 - data1[!,:zeta]
rel_error = abs.( error1 ) ./ abs.( data1[!,:zeta] )
maximum(rel_error)

