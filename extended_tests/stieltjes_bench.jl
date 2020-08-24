# check the generlised Euler (Stieltjes) numbers
#    this is mostly just a check that code is correct
#    most errors are arising because of errors in the functions being used
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot

S = Symbol("S_n")
data1 = CSV.read("../data/stieltjes_test_data_1.csv")
m = size(data1,1)
S1 = zeros(Float64, m)
S1_ref = zeros(Float64, m)
error1 = zeros(Float64, m)
rel_error1 = zeros(Float64, m)
for i=1:m
    S1[i] = stieltjes( data1[i,:n] )
    S1_ref[i] = data1[i,S]
    error1[i] = Float64( S1[i] - S1_ref[i]  )
    rel_error1[i] = error1[i]/Float64( S1_ref[i] )
end
println("max abs. error      = $(maximum( abs.(error1) ))")
println("max rel. abs. error = $(maximum( abs.(rel_error1) ))")
