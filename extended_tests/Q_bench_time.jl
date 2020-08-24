# compare times of direct calculation of c's compared to Crandall's recursive calculation
using Polylogarithms
using BenchmarkTools
using DataFrames, CSV
using PyPlot
using SpecialFunctions
import Base: parse
using Printf

Q = Polylogarithms.Q
c = Polylogarithms.c
max_terms = 10
max_n = 5

s = 5.0001;
r = s - n;
μ = -0.6931471805599453 + 3.141592653589793im #  = log(z) = log(-0.5) 
ℒ = log(complex(-μ))  # '\u2112'
 
result_direct = Array{Any}(undef, max_n, max_terms)
result_crandall = Array{Any}(undef, max_n, max_terms)

for k=0:max_n-1
    for j=1:max_terms
       if j<=3
            # b_direct = @benchmarkable Polylogarithms.c($k, $j, $ℒ)
            # result_crandall[k+1, n_terms] = run(b_direct)
            result_crandall[k+1, n_terms] =  minimum(@benchmark Polylogarithms.c($k, $j, $ℒ))
        end
        # b_crandall = @benchmarkable Polylogarithms.c_crandall($k, $j, $ℒ)
        # result_crandall[k+1, n_terms] = run(b_crandall)
        result_crandall[k+1, n_terms] = minimum(@benchmark Polylogarithms.c_crandall($k, $j, $ℒ))
   end
end
