
using Plots; plotlyjs()
default(show = true) # shouldn't need display commands

using SpecialFunctions
const SF = SpecialFunctions

SF.harmonic(1) - SF.gen_euler(0) - SF.digamma(2)
SF.harmonic(2) - SF.gen_euler(0) - SF.digamma(3)
SF.harmonic(3) - SF.gen_euler(0) - SF.digamma(4)


w = -log(complex(1.5))
lmu = log(-w)
lmu2 = log(w)
step = 0.002

tmp = 10.^collect(-6.0: 0.1: -1.0)
tmp = 10.^collect(-6.0: 0.1: -3.0)
x = vcat( -flipdim(tmp,1), 0.0, tmp)
plt1 = plot( x, real( SF.zeta.(1+x) + SF.gamma(-x).* w.^x ) ,
            title="n=0")
            scatter!( [0.0], [real(-lmu)] )

n = 0

# from Wolfram alpha, for n=0, https://www.wolframalpha.com/input/?i=zeta(1%2Bx)+%2B+gamma(-x)*w%5Ex
Q1 = (- SF.gen_euler(1)
      - 0.5*lmu2^2
      - γ*lmu2
      - π^2/12
      - γ^2/2
      )
# digamma(1) = - gamma

Q0 = 0 -lmu2
scatter!( [0.0], [real(Q0)] )
Q1 = -SF.gen_euler(1) - (SF.digamma(n+1)-lmu2)^2/2 - pi.^2/6 + SF.polygamma(1,n+1)/2
# pure wood: Q1 = SF.gen_euler(1) - (SF.digamma(n)-lmu)^2/2 + pi.^2/6 + SF.polygamma(1,n)/2
# Q2 = SF.gen_euler(2) - (SF.digamma(n)-lmu)^3/6 + (SF.digamma(n)-lmu)*( (1+x).^2/6 + SF.trigamma(n)/2) +  SF.polygamma(2,n)/6
Q2 =  SF.gen_euler(2) + (SF.digamma(1)-lmu2)^3/6 + [pi^2/6 + SF.polygamma(1,1)/2]*(SF.digamma(1)-lmu2) + SF.polygamma(2,1)/6
correction1 = Q0 + x .* Q1
correction2 = Q0 + x .* Q1  + x.^2 .* Q2 
plot!(x,  real(correction1), color=:red, linestyle=:dash) 
plot!(x,  real(correction2), color=:green, linestyle=:dash) 


tmp = NaN*zeros(size(x))
max_n = 5
a0 = zeros(max_n)
a1 = zeros(max_n)
a2 = zeros(max_n)
for n=1:max_n
    for i=1:length(tmp)
        if x[i]!=0.0
            tmp[i] = SF.gamma(-n-x[i])
        end
    end
    plt1 = plot( x, real( SF.zeta.(1+x) + (-1)^(n) * tmp .* (w).^(x) * factorial(n) ),
                reuse=false,
                title="n=$n")
    Q0 = SF.harmonic(n)-lmu2
    Q0 = SF.gen_euler(0) + SF.digamma(n+1)-lmu2
    scatter!( [0.0], [real(Q0)] )
    # Wood: Q1 = SF.gen_euler(1) - (SF.digamma(n+1)-lmu)^2/2 + (n+1+x).^2/6 + SF.polygamma(1,n+1)/2
    Q1 = -SF.gen_euler(1) - (SF.digamma(n+1)-lmu2)^2/2 - pi^2/6 + SF.polygamma(1,n+1)/2
    Q1 = -SF.gen_euler(1) - SF.digamma(n+1)^2/2 + SF.digamma(n+1)*lmu2 - lmu2^2/2 - pi^2/6 + SF.polygamma(1,n+1)/2
    a0[n] = SF.digamma(n+1) 
    γ22 = γ^2/2
    pi212 = π^2/12
    a1[n] = γ22 + pi212
    a2[n] = -γ^3/6 - γ*pi212 + SF.polygamma(2,1)/6.0
    for i=1:n
        a1[n] += (SF.harmonic(i)- γ)/i
        a2[n] += ( (SF.harmonic(i)^2 + SF.harmonic(i,2))/2.0 - γ*SF.harmonic(i) + γ^2/2 + pi212)/i
    end
    Q0 = SF.gen_euler(0) + a0[n] - lmu2
    Q1 = -SF.gen_euler(1) - a1[n] + a0[n]*lmu2 - lmu2^2/2
    # woodi-ish Q2 =  SF.gen_euler(2) + (SF.digamma(n)-lmu2)^3/6 + (-pi^2/6 + SF.polygamma(1,n)/2)*(SF.digamma(n)-lmu2) + SF.polygamma(2,n)/6
    Q2 = SF.gen_euler(2) + a2[n] - a1[n]*lmu2 + a0[n]*lmu2^2/2 - lmu2^3/6
    # Q2 = -5.5
    correction1 = Q0 + x .* Q1
    correction2 = Q0 + x .* Q1  + x.^2 .* Q2
    plot!(x,  real(correction1), color=:red, linestyle=:dash) 
    plot!(x,  real(correction2), color=:green, linestyle=:dash) 
end







# function c(n, j, L)
#     return (-1)^j/gamma(j+1)   *   gamma(j)   - b(n, j+1, L)
# end

# function b(n, j, L)
#     total = 0.0
#     for p=0:j
#         for t=0:j
#             for q=0:j
#                 if p+t+q==1
#                     total += L^p/gamma(p+1) * gamma(t,1)/gamma(t+1) * (-1)^(t+q) * f(k,q)
#                 end
#             end
#         end
#     end
#     return total
# end

# function f(n,q)
#     if q==0
#         return 1
#     else
#         total = 0.0
#         for h=0:q
#             total += f(n-1, q-h)
#         end
#     end
# end
