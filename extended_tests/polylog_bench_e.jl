# Check polylog
#    Compare m-th root and asymptotic series for very large z
# 
using Polylogarithms
using SpecialFunctions
using DataFrames, CSV
using PyPlot
using Printf
include("../src/utilities.jl")

import Polylogarithms.polylog_root
import Polylogarithms.polylog_asympt_series_1
L = Symbol("Li_s(z)")
    
  
# input data from Mathematica and reparse into complex numbers
C = 1
filename = @sprintf("../data/polylog_test_data_e_%d.csv", C)
data1 = CSV.read(filename, DataFrame; delim=",", comment="#", types=String)

C = 2
filename = @sprintf("../data/polylog_test_data_e_%d.csv", C)
data2 = CSV.read(filename, DataFrame; delim=",", comment="#", types=String)

#    has trouble reading in numbers like "2." so read all into strings, and parse
data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
data1[!,:r] = parse.(Float64, data1[!,:r] )
data1[!,:theta] = parse.(Float64, data1[!,:theta] )
data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
data1[!,L] = parse.(Complex{Float64}, data1[!,L] )

data2[!,:s] = parse.(Complex{Float64}, data2[!,:s] )
data2[!,:r] = parse.(Float64, data2[!,:r] )
data2[!,:theta] = parse.(Float64, data2[!,:theta] )
data2[!,:z] = parse.(Complex{Float64}, data2[!,:z] )
data2[!,L] = parse.(Complex{Float64}, data2[!,L] )

# check that complex numbers are being read in correctly
z = data1[!,:z]
z_2 = data2[!,:z]
z2 = data1[!,:r] .* exp.( im* data1[!,:theta] )
error = data1[!,:z] .- z2
maximum_abs_error = maximum(abs.(error))
# [data1[1:10,:z] parse.(Complex{Float64}, data1[1:10,:z] ) z2[1:10] ]
# [data1[1:10,:z]  z2[1:10] ]
# println("max error in parsing is ", maximum_abs_error)
max_x = log10( maximum( abs.( z ) ) )

# figure(1)
# plot(real.(z), imag.(z), "o")
# axis("equal")

m = size(data1,1)
Li = data1[!,L]
Li_b = data2[!,L]
s = data1[!,:s]
s_2 = data2[!,:s]
su = unique(s)
if length(su) > 1
    error()
else
    su = su[1]
    println("")
    @printf("s = %f + %f i\n", real(su), imag(su))
end
su_2 = unique(s_2)
if length(su_2) > 1
    error()
else
    su_2 = su_2[1]
    println("")
    @printf("s_2 = %f + %f i\n", real(su_2), imag(su_2))
end


S0 = zeros(Complex{Float64}, m)
n0 = zeros(Int64, m)
m0 = zeros(Int64, m)
error0 = zeros(Float64, m)
rel_error0 = zeros(Float64, m)

S1 = zeros(Complex{Float64}, m)
n1 = zeros(Int64, m)
error1 = zeros(Float64, m)
rel_error1 = zeros(Float64, m)

S1_b = zeros(Complex{Float64}, m)
n1_b  = zeros(Int64, m)
error1_b  = zeros(Float64, m)
rel_error1_b  = zeros(Float64, m)

S2 = zeros(Complex{Float64}, m)
n2 = zeros(Int64, m)
error2 = zeros(Float64, m)
rel_error2 = zeros(Float64, m)

for i=1:m
    # print(".")
    result0 = Polylogarithms.polylog_root(s[i], z[i])
    S0[i] = result0[1]
    n0[i] = result0[2]
    m0[i] = Int( floor(result0[3]/10) )
    error0[i] = abs( S0[i] - Li[i]  )
    rel_error0[i] = error0[i]/abs( Li[i] )

    result1 = Polylogarithms.polylog_asympt_series_1(s[i], z[i])
    S1[i] = result1[1]
    n1[i] = result1[2]
    error1[i] = abs( S1[i] - Li[i]  )
    rel_error1[i] = error1[i]/abs( Li[i] )
    
    result1_b = Polylogarithms.polylog_asympt_series_1(s_2[i], z_2[i])
    S1_b[i] = result1_b[1]
    n1_b[i] = result1_b[2]
    error1_b[i] = abs( S1_b[i] - Li_b[i]  )
    rel_error1_b[i] = error1_b[i]/abs( Li_b[i] )
    
    # result2 = - log(z[i])^s[i] / gamma(s[i]+1) # s neq -1,-2,-3,...
    result2 = - exp( s[i] * log(log(z[i])) - loggamma(s[i]+1) )# s neq -1,-2,-3,...
    S2[i] = result2
    n2[i] = 1
    error2[i] = abs( S2[i] - Li[i]  )
    rel_error2[i] = error2[i]/abs( Li[i] )

    
end
data1[:, :root_errors] = rel_error0
data1[:, :root_value] = S0
bad_data = filter( row -> row.root_errors > 1.0e-12, data1 )

fig = figure(@sprintf("../data/polylog_bench_e.csv"), figsize=(10,8))
clf()
title( @sprintf("s = %f + %f i", real(su), imag(su)) )
d = 0.025
ms = 3

subplot(221)
loglog( data1[!,:r] .- d, rel_error1, "go"; markersize=ms, label="asymptotic")
plot( data1[!,:r] .+ d, rel_error2, "rd"; markersize=ms, label="limit")
plot( data1[!,:r] .+ d, rel_error0, "b+"; markersize=ms, label="m-th root")
z = 10.0 .^ collect(2.0 : 1: 20 )
plot( z, 250.0 ./ z; linestyle="--", label="250 / |z|")
xlabel("|z|")
ylabel("relative absolute error")
legend()
ylim( 10 .^ [-16.0, 0] )
xlim( 10 .^ [2, max_x] )
xticks( 10 .^ collect(2.0 : 3 : max_x) )

subplot(223)
semilogx( data1[!,:r] .- d, n1, "go"; markersize=ms)
plot( data1[!,:r] .+ d, n2, "rd"; markersize=ms)
plot( data1[!,:r], n0, "b+"; markersize=ms)
    # θ = unique(data1[!,:theta])
    # for i=1:length(θ)
    #     k = findall(data1[!,:theta] .== θ[i])
    #     plot(data1[k,:r] .+ d, n2[k], "r-"; linewidth=0.5)
    #     text(data1[k[end],:r]+2*d, n2[k[end]], "$(θ[i]/π)"; verticalalignment="center", fontsize=9)
    # end
    # text(0.92, 65, "arg(z)/π"; verticalalignment="center")
xlabel("|z|")
ylabel("number of terms")
    # plot([0,1.1], [1,1]*Polylogarithms.default_max_iterations, "--")
ylim([0,70]) 
xlim( 10 .^ [2, max_x] )
xticks( 10 .^ collect(2.0 : 3 : max_x) )

# subplot(222)
# semilogy( data1[!,:theta] ./ π .- d, rel_error1, "go"; label="Series 4", markersize=ms)
# semilogy( data1[!,:theta] ./ π .+ d, rel_error2, "rd"; label="Limit", markersize=ms)
# semilogy( data1[!,:theta] ./ π, rel_error0, "b+"; label="Limit", markersize=ms)
# xlabel("arg(z)/π")
# ylabel("relative absolute error")
# # plot([-0.05,1.05], [1,1]*Polylogarithms.default_accuracy)
# # legend()

subplot(222)
loglog( data1[!,:r], rel_error1, "go"; markersize=ms, label="s = $(s[1])")
loglog( data2[!,:r], rel_error1_b, "rx"; markersize=ms, label="s = $(s_2[1])")
xlabel("|z|")
ylabel("relative absolute error")
legend()
ylim( 10 .^ [-16.0, 0] )
xlim( 10 .^ [2, max_x] )
xticks( 10 .^ collect(2.0 : 3 : max_x) )


subplot(224)
semilogx( data1[!,:r] , m0, "b+"; markersize=ms)
z = 10.0 .^ collect(2.0 : 1: max_x )
α = 1.2
plot(z, log.(z)/(pi*sqrt(α^2-1)); linestyle="--")
xlabel("|z|")
ylabel("m")
# plot([-0.05,1.05], [1,1]*Polylogarithms.default_max_iterations, "--")
ylim([0,30])
yticks( 0:5:30 )
xlim( 10 .^ [2, max_x] )
xticks( 10 .^ collect(2.0:3:max_x) )

savefig(@sprintf("plots/polylog_bench_e.pdf"); bbox_inches="tight")

