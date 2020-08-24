# comparison times for the zeta function
using SpecialFunctions
m = 100000
b = 1000.0
b = 8.0
b = 1.0

x = 2.0*b .* rand(m) .- b
y = 2.0*b .* rand(m) .- b
z = x + im*y
Z = Array{Complex{Float64}}(undef,m)

zeta(1.0) 
val, t, bytes, gctime, memallocs = @timed begin
    for i=1:m
        Z[i] = zeta(z[i])
    end
end
@printf("  finished %d calculations in time t=%.2f seconds, i.e., %.1f microsec each\n", m, t, (t/m)*1.0e6)



# Mathematica timing results
# z in [-1,1], [-8,8], [-16,16]
time_s = [16.06, 17.9, 18.9]
time_micros = 1.0e6 .* (time_s ./ 10000) # 



