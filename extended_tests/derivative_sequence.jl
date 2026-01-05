# check overall performance using randomly chosen points
# using Polylogarithms
# using DataFrames, CSV
using Plots
pyplot()

s = -4.0
s = -2.0 + 1.0*im
s = -2.0 + 0.1*im
s = -2.0

z = 0.5
z = 0.5 + 0.1 * im

max_k = 15
k = 1:max_k

sequence = z.^k  .* log.(k)  ./ k.^s

tester = log(z) .* k .* log.(k)    .-   s .* log.(k)  .+  1.0
ell = findall( abs.(tester) .> -1)

plot1 = plot(k, abs.(sequence); label="sequence", xlabel="k", ylabel="sequence value")
plot!( plot1, k[ell], abs.(tester[ell]); label="test" )









