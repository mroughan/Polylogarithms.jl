using SpecialFunctions
using Polylogarithms

"""
    zeta_derivative()

 created: 	Sun Jan  4 2026 
 email:   	matt@kanga2
 (c) Matt Roughan, 2026

 Calculates first derivative of the Riemann zeta function

## Arguments
* `n::Integer`: the number of elements to compute.
* `dim::Integer=1`: the dimensions along which to perform the computation.

## LateX
``\\alpha = \\beta``

## Examples
```jldoctest
julia> a = [1 2; 3 4]
2×2 Array{Int64,2}:
 1  2
 3  4
```
"""
function zeta_derivative( s::Number;
                          δ::Float64 = 0.2, 
                          )
    # convert to a floating point 64-bit type
    if isreal(s)
        s = convert(Float64, s)
    else
        s = convert(Complex{Float64}, s)
    end

    # choose the right series
    if abs(s - 1.0) < 1.0e-14
        return -Inf # zeta has simple pole at s=1
    elseif abs(s - 1.0) <= δ
        result = zeta_series_2( s; )

    elseif real(s) > 1 + δ
        result = zeta_series_1( s; )
        
    elseif real(s) < 0.0 - δ
        result = zeta_reflection( s; )
        
    elseif δ <= real(s) <= 1.0 + δ
        result = zeta_series_3( s; )
       
    else # real(s) near 0.0
        throw(DomainError(s, "haven't done this case yet"))
    end

    return result
end

max_log_ratio = 1001
log_ratios = log.(1:max_log_ratio) ./ log.(0:max_log_ratio-1)
    
# straight-forward calculation
#    zeta'(s) = - \sum_{k=2}^{\infty} ln(k) / k^s
function zeta_series_1( s;
                        max_iterations::Integer=10000000,
                        accuracy::Float64=1.0e-12,
                        )

    if real(s) <= 1
        throw(DomainError(s, "We should have real(s) > 1"))
    end
     
    total = 0.0
    converged = false
    k = 2
    a = -log(k) / k^s
    
    while k<=max_iterations && ~converged
        # println( " k = $k, a = $a ")
        total += a

        k = k+1
        a_old = a 
        a *= ( (k-1)/(k) )^s * log(k) ./ log(k-1)
        # a = -log(k) / k^s

        
        if abs(a)/abs(total) < 0.25*accuracy 
            converged = true
        end 
    end
    
    return total
end

# use this one near s=1.0
#    https://mathworld.wolfram.com/RiemannZetaFunction.html
#    https://math.stackexchange.com/questions/2513063/differential-of-zetas
# Near s=1
#    zeta(s) = 1/(s-1) + \sum_{n=0}^{\infty} (-1)^n γ_n (s-1)^n / n!
# then differentiate
# 
function zeta_series_2( s;
                        max_iterations::Integer=78, # max number of stieltjes constants available
                        accuracy::Float64=1.0e-12,
                        )

    total = -1 / (1-s)^2
    converged = false
    k = 1
    a = - 1.0
    b = stieltjes(1) * a
    ss = (s-1)
    
    while k<=max_iterations && ~converged
        total += b

        k = k+1
        a_old = a
        a *= ss / (k-1)
        b = a * stieltjes(k)
            
        if abs(b)/abs(total) < 0.25*accuracy 
            converged = true
        end 
    end
   
    return total
end

# use the functional equation
function zeta_reflection( s;
                          )

    if real(s) >= 0.0
        throw(DomainError(s, "We should have real(s) < 0"))
        # the functional equation is good for real(s) < 1, but we end up in no-mans land
    end
     
    term1 = log(2)
    term2 = log(π) 
    term3 = (π/2) * cot( s * π / 2)
    term4 = - polygamma(0, 1-s)
    term5 = zeta_derivative( 1 - s ) / zeta( 1- s)
    
    return zeta(s) * ( term1 + term2 + term3 + term4 + term5 )
end

# use Dirichlet eta function (see https://www.cecm.sfu.ca/~pborwein/PAPERS/P155.pdf)
#    zeta(s) = (1/ (1 - 2^(1-s)) ) \sum_{k=1}^{\infty} (-1)^(k+1) / k^s
# and differentiate
function zeta_series_3( s;
                        max_iterations::Integer=1000,
                        accuracy::Float64=1.0e-12,
                        )

    if real(s) <= 0.0
        throw(DomainError(s, "We should have real(s) > 0.0"))
    end

    total = 0.0
    converged = false
    k = 2
    a = log(k) / k^s
    
    while k<=max_iterations && ~converged
        total += a

        k = k+1
        a_old = a
        a *= - ( (k-1)/(k) )^s * log(k) ./ log(k-1)
        
        if abs(a)/abs(total) < 0.25*accuracy 
            converged = true
        end 
    end
    
    A =  1 + 2^(1-s) * log(2)
    B =  1 - 2^(1-s) 
    return total/B  + (B/A) * zeta(s)
end


# Riemann-Siegel algorithm for the critical line
#   http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
# but then need to differentiate this :( 

# there is also a Dirchlet series instead of the form used above
#   https://en.wikipedia.org/wiki/Riemann_zeta_function


# alternating series approach from 
#   http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
# still need to differentiate this
function zeta_series_4( s;
                        max_iterations::Integer=1000,
                        accuracy::Float64=1.0e-12,
                        )

    d = 12 # digits of 
    n = 1.3 * d + 0.9*imag(s) # so this is suitable for moderate imaginary part t

    d_n = zeros( n + 1 )
    for k = 0 : n

        a =
            
        for j = 0 : k-1
            d_n[j+1] += a
        end

    end
    d_n = n .* d_n
    
    B =  1 - 2^(1-s) 
    return total / (d_n[1] * B)
end
