# try out some ideas for sequences
#   starting with the direct sequence
using Polylogarithms
using DataFrames, CSV
using PyPlot
using Printf
include("utilities.jl")
using SpecialFunctions
using Base.MathConstants
twoπ = 2.0*π

# chose a "hard value"
s = complex(1.5)
z = 10.0 + 10.0im
x = log(-z)

result1 = polylog(s, z; min_iterations=1000 )
result2 = Polylogarithms.polylog_series_4( s, z )

max_n = 10
results = zeros(Complex{Float64}, max_n)
rel_errors = zeros(Float64, max_n)
for k = 0:1:max_n-1
    results[k+1] = (-1.0)^k * (1.0 - 2.0^(1-2*k)) * twoπ^(2*k) * bernoulli(2*k,0.0) * x^(s-2*k) / ( SpecialFunctions.gamma(2*k+1) * SpecialFunctions.gamma(s+1-2*k) )
    # rel_errors[k] = abs( results[k] - Li) / abs(Li)
end

a = diff(results)
a = results

figure()
k = 0:1:max_n-1
semilogy(k, abs.(a), "-")
plot(k, abs.(real.(a)), "-.")
plot(k, abs.(imag.(a)),"--")
# plot(k_star, abs.(a[k_star]), "o")



