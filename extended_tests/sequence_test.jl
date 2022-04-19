# try out some ideas for sequences
#   starting with the direct sequence
using Polylogarithms
using DataFrames, CSV
using PyPlot
using Printf
include(joinpath(@__DIR__, "..", "src", "utilities.jl"))
using SpecialFunctions

# chose a "hard value"
s = -7.46526-7.56069im
z1 = -7.60722-7.41506im
z = 1/z1
Li = 0.5783+114012.0im
result1 = polylog(s, z1; min_iterations=1000 )
result2 = Polylogarithms.polylog_series_2(s, z1; min_iterations=1000 )
error1 = Li - result1[1]
error2 = Li - result2[1]
rel_error1 = abs(error1) / abs(Li)
rel_error2 = abs(error2) / abs(Li)


k = 1:30
a = z.^k ./ k.^s
k_star = Int(ceil( real(s)/real(log(z)) ))

figure()
semilogy(k,abs.(a), "-")
plot(k, abs.(real.(a)), "-.")
plot(k, abs.(imag.(a)),"--")
plot(k_star, abs.(a[k_star]), "o")

let k=0
    tmp = z
    converged = false
    b = zeros(Complex{Float64}, size(a))
    max_iterations = length(a)
    min_iterations = k_star
    accuracy= 1.0e-12
    while k<max_iterations && ~converged
        k = k+1
        b[k+1] = b[k] + tmp
        tmp *= z * ( k/(k+1.0) )^s
        # println("   total = $total")
        if k > min_iterations && abs(tmp)/abs(b[k+1]) < 0.5*accuracy
            converged = true
        end
    end
    figure()
    plot(abs.(b[1:k]))
    println(b)
end

twoπ = 2*pi
G = (twoπ*im)^s * SpecialFunctions.zeta( 1-s, 0.5 + log(complex(-z1))/(twoπ*im) ) /  SpecialFunctions.gamma(s)
