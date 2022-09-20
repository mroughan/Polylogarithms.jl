# check overall performance using randomly chosen points
using Polylogarithms
using DataFrames, CSV
using PyPlot
using PyCall
using Printf
using StatsBase
fs = (4,4)
fs2 = (5,4)
PyCall.PyDict(matplotlib["rcParams"])["legend.markerscale"] = 3.0
# https://matplotlib.org/3.1.0/api/legend_api.html


L = Symbol("Li_s(z)")

# input data from Mathematica and reparse into complex numbers
C = 2
filename = @sprintf("../data/polylog_test_data_rand_%d.csv", C)
data1 = CSV.read(filename; delim=",", type=String)

# has trouble reading in numbers like "2." so read all into strings, and parse
data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
data1[!, L] = parse.(Complex{Float64}, data1[!, L] )

m = size(data1,1)
Li = data1[!, L]
s  = data1[!,:s]
z  = data1[!,:z]

# plot the locations of points
markerSize = 1

figure("z", figsize=fs)
title("distribution of z")
plot(real.(z), imag.(z), "."; markersize=markerSize, color=[0.8, 0.8, 0.8])
ylabel("Im(z)")
xlabel("Re(z)")
axis("equal")
    
figure("s", figsize=fs)
# title("distribution of s")
plot(real.(s), imag.(s), "."; markersize=markerSize, color=[0.8, 0.8, 0.8])
ylabel("Im(s)")
xlabel("Re(s)")
xticks(-8:2:8)
axis("equal")

S1 = zeros(Complex{Float64}, m)
S2 = zeros(Complex{Float64}, m)
n1 = zeros(Int64, m)
series1 = zeros(Int64, m)
error1 = zeros(Float64, m)
rel_error1 = zeros(Float64, m)

accuracy = 1.0e-12 # this is the default anyway, but can change it for testing
# polylog(1.5, 0.4)
polylog( complex(1.5), complex(0.4), Diagnostics() )
val, t, bytes, gctime, memallocs = @timed begin
    for i=1:m
        # print(".")
        S1[i],  n1[i], series1[i] = polylog(s[i], z[i], Diagnostics(); min_iterations=0, accuracy=accuracy )
        # S1[i] = result1[1]
        # n1[i] = result1[2]
        # series1[i] = result1[3]
    end
end
@printf("  finished %d calculations in time t=%.2f seconds, i.e., %.1f microsec each\n", m, t, (t/m)*1.0e6)

error1 = abs.( S1 - Li  )
rel_error1 = error1 ./ abs.( Li )
println("   max abs. rel. error1 = $(maximum( abs.(rel_error1) ))")

fig = figure(@sprintf("../data/polylog_bench_rand_1.csv"), figsize=(6,4))
clf()
h = hist(log10.(rel_error1), collect(-16 : 0.5 : -8); align="mid" )
ylabel("number")
xlabel("log10(relative absolute error)")
plot( log10(maximum(rel_error1))*[1,1], [0,maximum(h[1])], "--")
savefig(@sprintf("plots/polylog_bench_rand_%d.pdf", C); bbox_inches="tight")

figure("series", figsize=fs2)
ms = maximum(series1)
x = hist( series1, collect(0.5:1:ms+0.5); align="mid" )
for i=1:length(x[1])
    if i<=3
        @printf("  series %d used %d times\n", i, x[1][i])
    elseif i<8
        @printf("  reciprocal + series %d used %d times\n", i-3, x[1][i])
    elseif i< 20
        @printf("  duplication + series %d used %d times\n", i-10, x[1][i])
    else
        @printf("  duplication used more than once + series %d used %d times\n", i-10, x[1][i])
    end
end

# look at bad points
k = findall( rel_error1 .> 1.0e-12 )
count_bad = length(k)
bad_guys = DataFrame(s=s[k], z=z[k], Li=Li[k], err=rel_error1[k], series=series1[k],  n_terms=n1[k])

figure("z")
plot(real.(z[k]), imag.(z[k]), "r."; markersize=markerSize, color=[0.8, 0.0, 0.0])

figure("s")
plot(real.(s[k]), imag.(s[k]), "r."; markersize=markerSize, color=[0.8, 0.0, 0.0])


s_min = floor( minimum(min.(real.(s), imag.(s)) ))
s_max = ceil(  maximum(max.(real.(s), imag.(s)) ))
s_bins = collect(  s_min : 0.5 : s_max )
m = length(s_bins)
s_b = zeros(m-1)
error_v_s = zeros(m-1) 
for i=1:m-1
    s_b[i] = (s_bins[i] + s_bins[i+1])/2
    k_i = findall( s_bins[i] .<= imag.(s) .< s_bins[i+1] )
    error_v_s[i] = log10.( mean( rel_error1[k_i] ) )
end
figure("error vs s", figsize=fs2)
plot( s_b, error_v_s)
ylabel("number")
xlabel("series (values >4 are a sum of two components)")

# look at really bad points
k = findall( (rel_error1 .> 1.0e-8) .| isinf.(S1) )
count_really_bad = length(k)
really_bad_guys = DataFrame(s=s[k], z=z[k], Li=Li[k], err=rel_error1[k], series=series1[k],  n_terms=n1[k])

figure("z")
plot(real.(z[k]), imag.(z[k]), "r."; markersize=markerSize+3)
savefig( @sprintf("plots/polylog_bench_rand_errors_z_%d.pdf", C); bbox_inches="tight")

figure("s")
plot(real.(s[k]), imag.(s[k]), "r."; markersize=markerSize+3)
savefig( @sprintf("plots/polylog_bench_rand_errors_s_%d.pdf", C); bbox_inches="tight")


# show which series is used where
figure("z2", figsize=fs)
for i=1:ms
    k = findall( series1 .== i )
    if length(k) > 0
        if i<10
            plot(real.(z[k]), imag.(z[k]), "."; markersize=markerSize, label=@sprintf("Series %d",i))
        elseif i<50
            n_series = (i - floor(i,digits=-1)) / 2 # not right when we get to 50
            plot(real.(z[k]), imag.(z[k]), "."; markersize=markerSize, label=@sprintf("Series 2 used %d times", n_series))
        elseif i>=50
            i -= 50
            n_series = 5
            n_series += (i - floor(i,digits=-1)) / 2 # not right when we get to 50
            plot(real.(z[k]), imag.(z[k]), "."; markersize=markerSize, label=@sprintf("Series 2 used %d times", n_series))
        end
    end
end
legend()
ylabel("Im(z)")
xlabel("Re(z)")
axis("equal")
savefig( @sprintf("plots/polylog_bench_rand_domains_%d.pdf", C); bbox_inches="tight")

