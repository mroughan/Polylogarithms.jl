# check polylog
#    compare Series~2 and Series~3 near s=n > 0
#
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot
using Printf
include("utilities.jl")

series_1 = Polylogarithms.polylog_series_1
series_2 = Polylogarithms.polylog_series_2
series_3 = Polylogarithms.polylog_series_3
L = Symbol("Li_s(z)")
    
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
data1 = CSV.read("../data/polylog_test_data_b_1.csv"; delim=",", type=String)
#    has trouble reading in numbers like "2." so read all into strings, and parse
data1[!,:n] = Int.(parse.(Float64, data1[!,:n] ))
data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
data1[!,:r] = parse.(Float64, data1[!,:r] )
data1[!,:theta] = parse.(Float64, data1[!,:theta] )
data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
data1[!,L] = parse.(Complex{Float64}, data1[!,L] )
# check that complex numbers are being read in correctly
z = data1[!,:z]
z2 = data1[!,:r] .* exp.( im* data1[!,:theta] )
error = data1[!,:z] .- z2
maximum_abs_error = maximum(abs.(error))
# [data1[1:10,:z] parse.(Complex{Float64}, data1[1:10,:z] ) z2[1:10] ]
# [data1[1:10,:z]  z2[1:10] ]
println("max error in parsing is ", maximum_abs_error)
θ = unique(data1[!,:theta])
r = unique(data1[!,:r])

# figure(1)
# plot(real.(z), imag.(z), "o")
# axis("equal")

m = size(data1,1)
Li = data1[!,L]
s = data1[!,:s]

S1 = zeros(Complex{Float64}, m)
n1 = zeros(Int64, m)
error1 = zeros(Float64, m)
rel_error1 = zeros(Float64, m)

S2 = zeros(Complex{Float64}, m)
n2 = zeros(Int64, m)
error2 = zeros(Float64, m)
rel_error2 = zeros(Float64, m)

S3 = zeros(Complex{Float64}, m)
n3 = zeros(Int64, m)
error3 = zeros(Float64, m)
rel_error3 = zeros(Float64, m)

S5 = zeros(Complex{Float64}, m)
n5 = zeros(Int64, m)
error5 = zeros(Float64, m)
rel_error5 = zeros(Float64, m)

benchmark = series_3(data1[1,:n], z[1]; n_terms=3)[1]
bench_error = abs( benchmark - Li[1] )
bench_rel_error = bench_error / abs( Li[1] )
for i=1:m
    print(".")
    # result1 = series_1(s[i], z[i])
    # S1[i] = result1[1]
    # n1[i] = result1[2]
    # error1[i] = abs( S1[i] - Li[i]  )
    # rel_error1[i] = error1[i]/abs( Li[i] )

    result2 = series_2(s[i], z[i])
    S2[i] = result2[1]
    n2[i] = result2[2]
    error2[i] = abs( S2[i] - Li[i]  )
    rel_error2[i] = error2[i]/abs( Li[i] )

    result3 = series_3(s[i], z[i]; n_terms=3)
    S3[i] = result3[1]
    n3[i] = result3[2]
    error3[i] = abs( S3[i] - Li[i]  )
    rel_error3[i] = error3[i]/abs( Li[i] )

    result5 = series_3(s[i], z[i]; n_terms=5)
    S5[i] = result5[1]
    n5[i] = result5[2]
    error5[i] = abs( S5[i] - Li[i]  )
    rel_error5[i] = error5[i]/abs( Li[i] )
end
println("")
# println("max abs. error1 = $(maximum( abs.(error1) ))")
println("max abs. error2 = $(maximum( abs.(error2) ))")
println("max abs. error3 = $(maximum( abs.(error3) ))")
println("max abs. error5 = $(maximum( abs.(error5) ))")

fig = figure("polylog_bench_b", figsize=(6,4))
clf()
d = 0.00
 
# loglog( data1[!,:r] .- d, rel_error1, "o")
loglog( data1[!,:r] .+ d, rel_error2, "+"; label="Series 2" )
plot( data1[!,:r] .+ d, rel_error5, "rd"; label="Series 3 (5 terms)" )
plot( data1[!,:r] .+ d, rel_error3, "gx"; label="Series 3 (3 terms)" )
# plot( [minimum(r), maximum(r)], [1,1].*bench_rel_error, "--")
legend()
xlabel("|τ|")
ylabel("relative absolute error")
plot([minimum(r),maximum(r)], [1,1]*Polylogarithms.default_accuracy)
# xlim([0, 1.15])
# ylim([1.0-e15, 1.0e-10])
savefig("plots/polylog_bench_b.svg")
savefig("plots/polylog_bench_b.pdf")

# subplot(222)
# # loglog( data1[!,:r] .- d, rel_error1, "o")
# semilogy( data1[!,:r] .+ d, rel_error2, "+")
# plot( data1[!,:r] .+ d, rel_error3, "d")
# plot( [minimum(r), maximum(r)], [1,1].*bench_rel_error, "--")
# xlabel("r")
# ylabel("relative absolute error")
# plot([minimum(r),maximum(r)], [1,1]*Polylogarithms.default_accuracy)
# # xlim([0, 1.15])
# # ylim([1.0-e15, 1.0e-10])



# subplot(223)
# # loglog( data1[!,:r] .- d, n1, "o")
# loglog( data1[!,:r] .+ d, n2, "+")
# loglog( data1[!,:r] .+ d, n3, "d")
# # for i=1:length(θ)
# #     k = findall(data1[!,:theta] .== θ[i])
# #     plot(data1[k,:r] .+ d, n2[k], "g:")
# #     text(data1[k[end],:r]+2*d, n2[k[end]], "$(θ[i]/π)"; verticalalignment="center")
# # end
# # text(1.0, 72, "θ/π"; verticalalignment="center")
# xlabel("r")
# ylabel("number of terms")
# plot([minimum(r),maximum(r)], [1,1]*Polylogarithms.default_max_iterations)
# # xlim([0, 1.15])

# subplot(222)
# loglog( data1[!,:theta] ./ π .- d, rel_error1, "o")
# loglog( data1[!,:theta] ./ π .+ d, rel_error2, "+")
# xlabel("θ/π")
# ylabel("relative absolute error")
# plot([0,1], [1,1]*Polylogarithms.default_accuracy)

# subplot(224)
# loglog( data1[!,:theta] ./ π .- d, n1, "o")
# loglog( data1[!,:theta] ./ π .+ d, n2, "+")
# xlabel("θ/π")
# ylabel("number of terms")
# plot([0,1], [1,1]*Polylogarithms.default_max_iterations)
