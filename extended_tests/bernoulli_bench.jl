# check the Bernoulli numbers and polynomials
#    this is mostly just a check that code is correct
#    most errors are arising because of errors in the functions being used
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot
include("utilities.jl")

##############################################
# Bernoulli numbers
B = Symbol("B_n")
data1 = CSV.read("../data/bernoulli_test_data_1.csv", DataFrame)
m = size(data1,1)
B1 = zeros(Rational{Int64}, m)
B1_ref = zeros(Rational{Int64}, m)
error1 = zeros(Float64, m)
rel_error1 = zeros(Float64, m)
for i=1:m
    B1[i] = bernoulli( data1[i,:n] )
    B1_ref[i] = parse(Rational{Int64},data1[i,B])
    error1[i] = Float64( B1[i] - B1_ref[i]  )
    rel_error1[i] = error1[i]/Float64( B1_ref[i] )
end
println("max abs. error in Bernoulli numbers = $(maximum( abs.(error1) ))")

##############################################
# Bernoulli polynomials: identities
m = 20
B0 = zeros(Float64, m)
B0_ref = zeros(Rational{Int64}, m)
error0 = zeros(Float64, m)
for i=1:m
    B0_ref[i] = bernoulli( i )
    B0[i] = bernoulli( i, 0.0 )
    error0[i] = Float64( B0[i] - B0_ref[i]  )
end
figure(0)
clf()
semilogy( 1:m, abs.(error0) )
xlabel("n, (x=0.0)")
xticks(1:m)
ylabel("absolute error")
# ylim( 1.0e-17, 1.0e-15)
title("Bernoulli polynomials")
savefig("plots/bernoulli_bench_0.pdf")

##############################################
# Bernoulli polynomials: data v x
B = Symbol("B_n(x)")
data2 = CSV.read("../data/bernoulli_test_data_2.csv", DataFrame; delim=",", type=String)
data2[!,:n] = parse.(Float64, data2[!,:n] )
data2[!,:x] = parse.(Float64, data2[!,:x] )
data2[!,B] = parse.(Float64, data2[!,B] )
m = size(data2,1)
B2 = zeros(Float64, m)
B2_ref = zeros(Float64, m)
error2 = zeros(Float64, m)
rel_error2 = zeros(Float64, m)
for i=1:m
    B2[i] = bernoulli( Int(data2[i,:n]),  data2[i,:x])
    B2_ref[i] = data2[i,B]
    error2[i] = Float64( B2[i] - B2_ref[i]  )
    rel_error2[i] = error2[i] / B2_ref[i]
end
figure(2)
clf()
semilogy( data2[!,:x], abs.(rel_error2) )
xlabel("x, (n=12)")
ylabel("absolute rel. error")
# ylim( 1.0e-17, 1.0e-15)
title("Bernoulli polynomials")
savefig("plots/bernoulli_bench_2.pdf")


##############################################u
# Bernoulli polynomials: data v n
B = Symbol("B_n(x)")
data3 = CSV.read("../data/bernoulli_test_data_3.csv", DataFrame)
m = size(data3,1)
B3 = zeros(Float64, m)
B3_ref = zeros(Float64, m)
error3 = zeros(Float64, m)
rel_error3 = zeros(Float64, m)
for i=1:m
    B3[i] = bernoulli( Int(data3[i,:n]),  data3[i,:x])
    B3_ref[i] = data3[i,B]
    error3[i] = Float64( B3[i] - B3_ref[i]  )
    rel_error3[i] = error3[i] / B3_ref[i]
end
figure(3)
clf()
semilogy( data3[!,:n], abs.(rel_error3) )
xlabel("n, (x=1.5)")
ylabel("absolute rel. error")
# ylim( 1.0e-17, 1.0e-15)
title("Bernoulli polynomials")
savefig("plots/bernoulli_bench_3.pdf")



