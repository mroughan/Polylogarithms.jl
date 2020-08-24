# check the Bernoulli numbers and polynomials
#    this is mostly just a check that code is correct
#    most errors are arising because of errors in the functions being used
using Polylogarithms
using PyPlot
using SpecialFunctions
import Base: parse
using Printf

Q = Polylogarithms.Q
c1 = Polylogarithms.c
c2 = Polylogarithms.c_crandall
b = Polylogarithms.b_crandall
f = Polylogarithms.f_crandall
g = Polylogarithms.g_crandall

# check f
@printf("Testing f\n")
for k=0:5
    @printf(" k = %d\n", k)
    @printf("   f(k,1) = %f, -H_k= %f\n", f(k,1), -harmonic(k) )
    @printf("   f(k,2) = %f, (H_k^2 + H_{k,2})/2= %f\n", f(k,2), (harmonic(k)^2 + harmonic(k,2))/2 )
end
@printf("\n")


# # check g against numerical derivates of Gamma at 1
@printf("Testing g\n")
function G(t::Integer, x::Float64; dx::Float64 = 1.0e-3)
    if t==0
        return 1
    elseif t==1
        return ( SpecialFunctions.gamma(x+dx) - SpecialFunctions.gamma(x) )/dx
    elseif t>1
        return ( G(t-1,x+dx) - G(t-1,x) ) /dx 
    else
        throw(DomainError(t))
    end
end
for t=0:5
    @printf(" t = %d, ", t)
    @printf(" G^{t}(1) = %.7f,  g(t) = %.7f\n",  G(t,1.0), g(t) ) 
end
@printf("\n")

# check c
z = complex(-1.0/2.0)
μ = log(z)
ℒ = log(-μ)
for k=0:5
    @printf(" k = %d\n", k)
    for j=0:2
        @printf("  j = %d,  ", j)
        @printf("  c_explicit = %.12f + %.12fi, ", real(c1(k, j, ℒ)), imag(c1(k, j, ℒ)) )
        @printf("  c_crandall = %.12f + %.12fi, ", real(c2(k, j, ℒ)), imag(c2(k, j, ℒ)) )
        @printf("\n")
   end
   @printf("\n")
end

c2(1, 3, ℒ)
