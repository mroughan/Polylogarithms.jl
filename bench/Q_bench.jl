# check the Bernoulli numbers and polynomials
#    this is mostly just a check that code is correct
#    most errors are arising because of errors in the functions being used
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot
using SpecialFunctions
import Base: parse
using Printf

Q = Polylogarithms.Q
c = Polylogarithms.c
# b = Polylogarithms.b
# f = Polylogarithms.f
# g = Polylogarithms.g
    
# from https://github.com/JuliaLang/julia/issues/18328
function parse(::Type{Complex{T}}, s::AbstractString) where {T<:Real}
    s = replace(s, "*^" => "e") # Mathematica seems to export in forms like "3.061616997868383*^-18"
    real_r = r"([+-]?\d*\.\d*(e[+-]\d+)?)([\s+-]|$)"
    imag_r = r"([+-]?\d*\.\d*(e[+-]\d+)?)\s*\*\s*I"
    
    z = zero(Complex{T})
    
    m = match( imag_r, s )
    if m != nothing
        z += im*parse(T, m[1])
    end
    
    m = match( real_r, s )
    if m != nothing
        z += parse(T, m[1])
    end
    return z
end

# input data from Mathematica and reparse into complex numbers
data1 = CSV.read("../data/Q_test_data_1.csv"; delim=",", type=String)
#    has trouble reading in numbers like "2." so read all into strings, and parse
data1[!,:n] = Int.(parse.(Float64, data1[!,:n] ))
data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
data1[!,:r] = parse.(Float64, data1[!,:r] )
data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
data1[!,:Q] = parse.(Complex{Float64}, data1[!,:Q] )
n = data1[!,:n]
m = size(data1,1)
s = data1[!,:s]
r = data1[!,:r]
μ = log.( data1[!,:z] )
ℒ = log.(complex.(-μ))  # '\u2112'

# figure(1)
# semilogx( data1[!,:r], real( data1[!,:Q] ), "b" )
# semilogx( data1[!,:r], imag( data1[!,:Q] ), "r--" )

# figure(2)
# loglog( data1[!,:r], abs.(real( data1[!,:Q] .- data1[1,:Q] )), "b" )
# loglog( data1[!,:r], abs.(imag( data1[!,:Q] .- data1[1,:Q] )), "r--" )

n_terms = 7

result = zeros(Complex{Float64}, m, n_terms)
error = zeros(Float64, m, n_terms)
rel_error = zeros(Float64, m, n_terms)

result2 = zeros(Complex{Float64}, m)
error2 = zeros(Float64, m)
rel_error2 = zeros(Float64, m)

for i=1:m
    for j=1:n_terms
        print("$i.")
        result[i,j] = Q(n[i], r[i], ℒ[i]; n_terms = j)
        error[i,j] = abs( result[i,j] - data1[i,:Q]  )
        rel_error[i,j] = error[i,j]/abs( data1[i,:Q] )
         
        result2[i] = (-1)^n[i] * factorial(n[i]) * gamma(-n[i] .- r[i]) * (-μ[i])^r[i]  + zeta(1 + r[i])
        error2[i] = abs( result2[i] - data1[i,:Q]  )
        rel_error2[i] = error2[i]/abs( data1[i,:Q] )
    end
end
println("")

for i=1:5
    @printf(" r1 = %.20f, r2 = %.20f, r3 = %.20f\n", abs(result[i,1]), abs(result[i,2]) , abs(result[i,3])  )
end


fig = figure("Q_bench", figsize=(5,5))
clf()
d = 0.00
for j=1:n_terms
    loglog( data1[!,:r] .+ d, rel_error[:,j], "-"; label=@sprintf("%d terms",j) )
end
loglog( data1[!,:r] .+ d, rel_error2[:], "r-"; label="direct")
legend()
xlabel("tau")
ylabel("relative absolute error")
