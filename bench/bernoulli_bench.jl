# check the Bernoulli numbers and polynomials
#    this is mostly just a check that code is correct
#    most errors are arising because of errors in the functions being used
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot
import Base: parse
# from https://github.com/JuliaLang/julia/issues/18328
function parse(::Type{Rational{T}}, x::AbstractString) where {T<:Integer}
    list = split(x, '/', keepempty=false)
    if length(list) == 1
        return parse(T, list[1]) // 1
    else
        @assert length(list) == 2
        return parse(T, list[1]) // parse(T, list[2])
    end
end

##############################################
# Bernoulli numbers
B = Symbol("B_n")
data1 = CSV.read("../data/bernoulli_test_data_1.csv")
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
println("max abs. error = $(maximum( abs.(error1) ))")


##############################################
# Bernoulli polynomials: identities
m = 20
B2 = zeros(Float64, m)
B2_ref = zeros(Rational{Int64}, m)
error2 = zeros(Float64, m)
for i=1:m
    B2_ref[i] = bernoulli( i )
    B2[i] = bernoulli( i, 0.0 )
    error2[i] = Float64( B2[i] - B2_ref[i]  )
end
figure(2)
clf()
semilogy( 1:m, abs.(error2) )
xlabel("n, (x=0.0)")
ylabel("absolute error")
# ylim( 1.0e-17, 1.0e-15)
title("Bernoulli polynomials")
savefig("plots/bernoulli_bench_2.pdf")

##############################################
# Bernoulli polynomials: data v x
B = Symbol("B_n(x)")
data3 = CSV.read("../data/bernoulli_test_data_2.csv")
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
semilogy( data3[!,:x], abs.(rel_error3) )
xlabel("x, (n=12)")
ylabel("absolute rel. error")
# ylim( 1.0e-17, 1.0e-15)
title("Bernoulli polynomials")
savefig("plots/bernoulli_bench_3.pdf")


##############################################
# Bernoulli polynomials: data v n
B = Symbol("B_n(x)")
data4 = CSV.read("../data/bernoulli_test_data_3.csv")
m = size(data4,1)
B4 = zeros(Float64, m)
B4_ref = zeros(Float64, m)
error4 = zeros(Float64, m)
rel_error4 = zeros(Float64, m)
for i=1:m
    B4[i] = bernoulli( Int(data4[i,:n]),  data4[i,:x])
    B4_ref[i] = data4[i,B]
    error4[i] = Float64( B4[i] - B4_ref[i]  )
    rel_error4[i] = error4[i] / B4_ref[i]
end
figure(4)
clf()
semilogy( data4[!,:n], abs.(rel_error4) )
xlabel("n, (x=1.5)")
ylabel("absolute rel. error")
# ylim( 1.0e-17, 1.0e-15)
title("Bernoulli polynomials")
savefig("plots/bernoulli_bench_4.pdf")



