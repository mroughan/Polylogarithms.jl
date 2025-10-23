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
accuracy = 1.0e-12 # this is the default anyway, but can change it for testing

# input data from Mathematica and reparse into complex numbers
error1 = Vector{Vector{Float64}}(undef,4)
rel_error1 = Vector{Vector{Float64}}(undef,4)
for C=1:4
    filename = @sprintf("../data/polylog_test_data_rand_%d.csv", C)
    data1 = CSV.read(filename, DataFrame; delim=",", types=String)
    
    # has trouble reading in numbers like "2." so read all into strings, and parse
    data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
    data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
    data1[!, L] = parse.(Complex{Float64}, data1[!, L] )
    
    m = size(data1,1)
    Li = data1[!, L]
    s  = data1[!,:s]
    z  = data1[!,:z]
#     error1[C] = Vector{Float64}(undef,m)
#     rel_error1[C] = Vector{Float64}(undef,m)
    
    S1 = zeros(Complex{Float64}, m)
    n1 = zeros(Int64, m)
    series1 = Vector{Any}(undef,m)
    max_recursion1 = zeros(Int64, m)

    # polylog(1.5, 0.4)
    polylog( complex(1.5), complex(0.4), Diagnostics() )
    val, t, bytes, gctime, memallocs = @timed begin
        for i=1:m
            # print(".")
            clearcache()
            S1[i],  n1[i], series1[i], max_recursion1[i] = polylog(s[i], z[i], Diagnostics(); min_iterations=0, accuracy=accuracy )
        end
    end
    println("C = $C")
    @printf("    finished %d calculations in time t=%.2f seconds, i.e., %.1f microsec each\n", m, t, (t/m)*1.0e6)
 
    error1[C] = abs.( S1 - Li  ) 
    rel_error1[C] = error1[C] ./ abs.( Li )
    println("     max abs. rel. error1 = $(maximum( rel_error1[C] ))")
    println("     no. of points outside bounds = ", length(findall( rel_error1[C] .>  accuracy) ))
end

fig = figure(@sprintf("../data/polylog_bench_rand.csv"), figsize=(6,4))
clf()
h = hist([log10.(rel_error1[1]), log10.(rel_error1[2]), log10.(rel_error1[3]), log10.(rel_error1[3])], collect(-16 : 0.5 : -10);
         align="mid", alpha=0.7, label=["dataset 1", "dataset 2", "dataset 3", "dataset 3"],
         )
# hist(log10.(rel_error1[2]), collect(-16 : 0.5 : -10); align="mid", alpha=0.3 )
# hist(log10.(rel_error1[3]), collect(-16 : 0.5 : -10); align="mid", alpha=0.3 )
ylabel("number")
xlabel("log10(relative absolute error)")
legend()
# plot( log10(maximum(rel_error1))*[1,1], [0,maximum(h[1])], "--")
savefig("plots/polylog_bench_rand.pdf"; bbox_inches="tight")
 
