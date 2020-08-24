# check the Harmonic numbers
#    this is mostly just a check that code is correct
#    most errors are arising because of errors in the functions being used
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot

H = Symbol("real(H_n)")

##############################################
# integer arguments
data1 = CSV.read("../data/harmonic_test_data_1.csv", DataFrame)
m = size(data1,1)
H1 = zeros(Float64, m)
error1 = zeros(Float64, m)
rel_error1 = zeros(Float64, m)
for i=1:m
    H1[i] = harmonic( data1[i,:n] )
    error1[i] = H1[i] -  data1[i,H]
    rel_error1[i] = error1[i]/data1[i,H]
end
figure(1)
clf()
semilogy( data1[!,:n], abs.(rel_error1) )
xlabel("n")
ylabel("relative absolute error")
ylim( 1.0e-17, 1.0e-15)
title("Harmonic numbers")
savefig("plots/harmonic_bench_1.pdf")

##############################################
# real arguments
data2 = CSV.read("../data/harmonic_test_data_2.csv", DataFrame)
m = size(data2,1)
H2 = zeros(Float64, m)
error2 = zeros(Float64, m)
rel_error2 = zeros(Float64, m)
for i=1:m
    H2[i] = harmonic( data2[i,:x] )
    error2[i] = H2[i] -  data2[i,H]
    rel_error2[i] = error2[i]/data2[i,H]
end
figure(2)
clf()
# plot( data2[!,:x], rel_error2 )
semilogy( data2[!,:x], abs.(rel_error2) )
k2 = findall( abs.(rel_error2) .> 1.0e-13 )
plot( data2[k2,:x], abs.(rel_error2[k2]), "r" )
xlabel("x")
ylabel("relative absolute error")
title("Harmonic numbers (real arguments)")
savefig("plots/harmonic_bench_2.pdf")

##############################################
# complex arguments
zr = Symbol("Re(z)")
zi = Symbol("Im(z)")
Hr = Symbol("Re(H(z))")
Hi = Symbol("Im(H(z))")
data3 = CSV.read("../data/harmonic_test_data_3.csv", DataFrame)
m = size(data3,1)
H3 = zeros(Complex{Float64}, m)
error3 = zeros(Complex{Float64}, m)
rel_error3 = zeros(Float64, m)
for i=1:m
    z = data3[i,zr] + im*data3[i,zi]
    Hz = data3[i,Hr] + im*data3[i,Hi]
    # println(" i = $i, z = $z, typeof(z)=$(typeof(z)), H = $(Hz)" )
    H3[i] = harmonic( z )
    error3[i] = H3[i] - Hz
    rel_error3[i] = abs(error3[i])/ abs( Hz )
end
figure(3)
clf()
# plot( data3[!,:x], rel_error3 )
semilogy( atan.(data3[!,zi],data3[!,zr])/pi, abs.(rel_error3) )
k3 = findall( abs.(rel_error3) .> 1.0e-13 )
plot( data3[k3,:x], abs.(rel_error3[k3]), "r" )
xlabel("θ/π") 
ylabel("relative absolute error")
title("Harmonic numbers (complex arguments)")
savefig("plots/harmonic_bench_3.pdf")


##############################################
# integer arguments for generalize H
H = Symbol("real(H_{n,r})")
data4 = CSV.read("../data/harmonic_test_data_4.csv", DataFrame)
m = size(data4,1)
H4 = zeros(Float64, m)
error4 = zeros(Float64, m)
rel_error4 = zeros(Float64, m)
for i=1:m
    H4[i] = harmonic( data4[i,:n], data4[i,:r] )
    error4[i] = H4[i] -  data4[i,H]
    rel_error4[i] = error4[i]/data4[i,H]
end
figure(4)
clf()
semilogy( data4[!,:n] + 0.1*data4[!,:r], abs.(rel_error4) )
xlabel("n + 0.1*r")
ylabel("relative absolute error")
ylim( 1.0e-17, 1.0e-15)
title("Generalised harmonic numbers")
savefig("plots/harmonic_bench_4.pdf")

##############################################
# integer n, real r for generalized H
H = Symbol("real(H_{n,r})")
data5 = CSV.read("../data/harmonic_test_data_5.csv", DataFrame)
m = size(data5,1)
H5 = zeros(Float64, m)
error5 = zeros(Float64, m)
rel_error5 = zeros(Float64, m)
for i=1:m
    H5[i] = harmonic( data5[i,:n], data5[i,:r] )
    error5[i] = H5[i] -  data5[i,H]
    rel_error5[i] = error5[i]/data5[i,H]
end
figure(5)
clf()
semilogy( data5[!,:n] + 0.1*data5[!,:r], abs.(rel_error5) )
xlabel("n + 0.1*r") 
ylabel("relative absolute error")
ylim( 1.0e-17, 1.0e-15)
title("Generalised harmonic numbers")
savefig("plots/harmonic_bench_5.pdf")
