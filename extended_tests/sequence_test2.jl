# try out some ideas for sequences
#   starting with the direct sequence
using Polylogarithms
using DataFrames, CSV
using PyPlot
using Printf
include(joinpath(@__DIR__, "..", "src", "utilities.jl"))
using SpecialFunctions

# chose a "hard value"
s = -6.084227492516975 - 7.500776841841313im
z = 1.4913747102091242 - 3.6642351854073443im
Li = 26.960114195028023 - 13.176237391440013im

μ = log(z)
(-μ)^(s-1) * SpecialFunctions.gamma(1-s)

result1 = polylog(s, z; min_iterations=1000 )
result2 = Polylogarithms.polylog_series_2(s, z; min_iterations=0 )
result2a = Polylogarithms.polylog_series_2(s, z; max_iterations=1 )
result3 = Polylogarithms.polylog_series_2(s, z; min_iterations=100 )
error1 = Li - result1[1]
error2 = Li - result2[1]
error3 = Li - result3[1]
rel_error1 = abs(error1) / abs(Li)
rel_error2 = abs(error2) / abs(Li)
rel_error3 = abs(error3) / abs(Li)

# a = z.^k ./ k.^s
# k_star = Int(ceil( real(s)/real(log(z)) ))

max_n = 100
results = zeros(Complex{Float64}, max_n)
rel_errors = zeros(Float64, max_n)
for k = 1:1:max_n
    results[k] = Polylogarithms.polylog_series_2(s, z; min_iterations=k, max_iterations=k)[1]
    rel_errors[k] = abs( results[k] - Li) / abs(Li)
end

a = diff(results)

figure()
k = 2:1:max_n
semilogy(k, abs.(a), "-")
plot(k, abs.(real.(a)), "-.")
plot(k, abs.(imag.(a)),"--")
# plot(k_star, abs.(a[k_star]), "o")



# let k=0
#     tmp = z
#     converged = false
#     b = zeros(Complex{Float64}, size(a))
#     max_iterations = length(a)
#     min_iterations = k_star
#     accuracy= 1.0e-12
#     while k<max_iterations && ~converged
#         k = k+1
#         b[k+1] = b[k] + tmp
#         tmp *= z * ( k/(k+1.0) )^s
#         # println("   total = $total")
#         if k > min_iterations && abs(tmp)/abs(b[k+1]) < 0.5*accuracy
#             converged = true
#         end
#     end
#     figure()
#     plot(abs.(b[1:k]))
#     println(b)
# end


