using SpecialFunctions
using Polylogarithms

const max_log_ratio = 1001
const log_ratios = log.(1:max_log_ratio) ./ log.(0:max_log_ratio-1)


"""
    zeta_derivative()

 Calculates first derivative of the Riemann zeta function

## Arguments
* `n::Integer`: the number of elements to compute.
* `δ::Float64 = 0.2`: the point at which we swap to the Laurent series around s=1.0
* `accuracy::Float64=1.0e-14`: intended accuracy, although this doesn't have as much affect here
* `n_terms::Integer=0`: number of terms in series 4 if used (default is zero lets it set this itself)

## Returns (a tuple)
* Result
* Number of terms used in series (or a diagnostic for non-series approaches)

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> zeta_derivative(0.5)
(-3.9226461392091503, 19)
```
"""
function zeta_derivative( s::Number;
                          δ::Float64 = 0.2, 
                          accuracy::Float64=1.0e-14,
                          n_terms::Integer=0,
                          )
    # convert to a floating point 64-bit type
    if isreal(s)
        s = convert(Float64, s)
    else
        s = convert(Complex{Float64}, s)
    end
   
    # choose the right approach
    if abs(s - 1.0) < 1.0e-14
        result = (-Inf, 0) # zeta has simple pole at s=1
    elseif abs(s - 1.0) <= δ
        result = zeta_series_2( s; accuracy=accuracy)
    elseif real(s) >= 1/2
        # alternating series approach
        result = zeta_series_4( s; accuracy=accuracy, n_terms=n_terms)
     elseif abs(s) <= δ
        # reflection doesn't work well near s=0 because of cancelling poles
        #   we should do our own Taylor series around 0, but the alternating series works OK near here
        result = zeta_series_4( s; accuracy=accuracy, n_terms=n_terms)
#    elseif n<0 && abs(s - 2*round( s/2 )) < 1.0e-6 
#       # old reflection doesn't work well near the trivial zeros at -2,-4, -6, ...
#       # but neither does the alternating series
#       #        result = zeta_d_near_zeros( s )
#       # but new reflection included both bits
    elseif real(s) < 1/2
        result = zeta_reflection2( s; δ=δ, accuracy=accuracy, n_terms=n_terms)
    else 
        # we really should do a separate case for large imag(s)
    end

    return result[1], result[2]
end
    
# straight-forward calculation
#    zeta'(s) = - \sum_{k=2}^{\infty} ln(k) / k^s
# which converges super slowly, so never use this
function zeta_series_1( s::Number;
                        max_iterations::Integer=10000000,
                        accuracy::Float64=1.0e-14,
                        )
    if real(s) <= 1
        throw(DomainError(s, "We should have real(s) > 1"))
    end
    
    total = 0.0
    converged = false
    k = 2
    a = -log(k) / k^s
    while k<=max_iterations && ~converged
        total += a
       k = k+1
        a_old = a 
        a *= ( (k-1)/(k) )^s * log(k) ./ log(k-1)
         if abs(a)/abs(total) < 0.25*accuracy 
            converged = true
        end 
    end
    return total, 1000000 + k
end

# Laurent series expansion of ζ(s) around s = 1
# use this one very near s=1.0
#    https://mathworld.wolfram.com/RiemannZetaFunction.html
#    https://math.stackexchange.com/questions/2513063/differential-of-zetas
#    https://mast.queensu.ca/~murty/ram-siddhi.pdf
# Near s=1
#    zeta(s) = 1/(s-1) + \sum_{n=0}^{\infty} (-1)^n γ_n (s-1)^n / n!
# then differentiate
# 
function zeta_series_2( s::Number;
                        max_iterations::Integer=78, # max number of stieltjes constants available
                        accuracy::Float64=1.0e-14,
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
    if k >= max_iterations
        error(" did not converge fast enough" )
    end
   
    return total, 1000 + k
end

# use the functional equation to calculate for s<1/2
#   COMPUTATION OF THE DERIVATIVES OF THE RIEMANN ZETA-FUNCTION IN THE COMPLEX DOMAIN, U.S.S.R. Comput.Maths.Math.Phys.,Vo1.28,No.4,pp.115-124,1988, A.YU. YEREMIN, I.E. KAPORIN and M.K. KERIMOV (equation 3.2)
#   WALTHER A., Anschauliches zur Riemannschen Zeta-Function.Acta Math., 46, 393-400, 1926.
#
# NOTE: this is the old one that doesn't work for negative (even) integers because of trivial zeros of ζ(-2m)
# 
function zeta_reflection( s::Number;
                         δ::Float64 = 0.2, 
                         accuracy::Float64=1.0e-14,
                         )

    if real(s) >= 0.5
        throw(DomainError(s, "We should have real(s) < 1/2"))
        # the functional equation is good for a awder range, but we shouldn't need it
    end

    term1 = log(2π)
    term2 = (π/2) * cot( s * π / 2)
    term3 = - polygamma(0, 1-s)
    result = zeta_derivative( 1 - s; δ=δ, accuracy=accuracy )
    term4 = - result[1] / zeta( 1 - s )
    
    return zeta(s) * ( term1 + term2 + term3 + term4 ), -result[2]
end


# incorporate
# COMPUTATION OF THE DERIVATIVES OF THE RIEMANN ZETA-FUNCTION IN THE COMPLEX DOMAIN, U.S.S.R. Comput.Maths.Math.Phys.,Vo1.28,No.4,pp.115-124,1988, A.YU. YEREMIN, I.E. KAPORIN and M.K. KERIMOV (equation 3.2)
#   equation 3.1 into 3.2 to avoid the poles
#
function zeta_reflection2( s::Number;
                           δ::Float64 = 0.2, 
                           accuracy::Float64=1.0e-14,
                           n_terms::Integer=0
                           )

    if real(s) >= 0.5
        throw(DomainError(s, "We should have real(s) < 1/2"))
        # the functional equation is good for a awder range, but we shouldn't need it
    end

    z1 = zeta(1-s)
    z2 = zeta(s)
    result = zeta_derivative( 1 - s; δ=δ, accuracy=accuracy, n_terms=n_terms )
    zd = result[1]

    g = gamma(1-s)
    pg = polygamma(0, 1-s)
    
    ct = (π/2) * cot( s * π / 2)
    cst = (π/2) * cos( s * π / 2)
    
    total = log(2*pi) * z2     +
               2.0 ^ (s-1) *  pi^s  * cos(pi*s/2) * g * z1    -
               pg * z2                                       -
               2.0 ^ s *  pi^(s-1)  * sin(pi*s/2) * g * zd 
        
    return total, -result[2]
end


# alternating series approach from 
#   http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
#   Series with binomial‑like coefficients for the Riemann zeta, Igoris Belova, s10231-021-01142-1.pdf
#  The Riemann Zeta-Function and Its Derivatives, Bejoy K. Choudhury, Choudhury-RiemannZetaFunctionDerivatives-1995.pdf
# uses Dirichlet eta function (see https://www.cecm.sfu.ca/~pborwein/PAPERS/P155.pdf)
#    zeta(s) = (1/ (1 - 2^(1-s)) ) \sum_{k=1}^{\infty} (-1)^(k+1) / k^s
# still need to differentiate this
function zeta_series_4( s::Number;
                       accuracy::Float64=1.0e-14,
                       n_terms::Integer=0, # number of terms
                       )

    d = ceil(-log10(accuracy))
       
    if n_terms==0
        # http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
        #     prolly needs adjustment when derivative is calculated
        n_terms = Int64( ceil(1.3*d + 0.9*abs(imag(s)) ) ) + 1 # added 1 here for big real(s) where we seem to have issues
    else
        # keep n
    end
    if n_terms > 300
        throw(DomainError(s, "|imag(s)| too big -- too many terms in sequences"))
    end
    
    m, d, u = d_calc_1( n_terms )
    
    total = 0.0 
    for k=1:n_terms
        total -= (-1)^(k-1) * m[k] * log(k) / k^s
    end

    de_etaing_term = 1 - 2^(1-s)
    de_etaing_derivative= -log(2) * 2^(1-s)
    result = (total + zeta(s)*de_etaing_derivative )/ de_etaing_term

    # note that n is chosen to set the absolute error, not relative
    # so once we have a result, check the relative error bound, and if neeeded recalc, but note that we can't do
    # much better than machine precision, so this won't work for very small values
    
    return result, n_terms
end

# I'm not using this at the moment, just the rule of thumb, but we might in the future
function zeta_error_bound( s::Number, n::Integer )
    # note this will behave a bit weird around s=1
    
    t = abs(imag(s))
    term1 = 3 / (3 + sqrt(8))^n 
    term2 = ( 1 + 2*t ) * exp( t * pi / 2)
    term3 = abs(  1 - 2^(1-s) )

    return term1 * term2 / term3
end

# direct and stupid calculation for testing small n
#   don't use this -- just keep it for comparisons
function d_calc_0( n::Integer )
    d = zeros(n+1) # note indexes offset, so d[1] = d_0

    k = 0:n
    u = n .* factorial.(n .+ k .- 1) .* (4 .^k) ./ ( factorial.(n .- k) .* factorial.(2 .* k) )
    
    d[1] = 1
    u[1] = 1
    for k=1:n
        d[k+1] = d[k] + u[k+1]
    end

    m = 1 .- d ./ d[n+1]

    return m, d, u
end

# better iterative approach
#   but we should precalculate d, or at least cache it as n will be the same as long as imag(s)=const
function d_calc_1( n::Integer )

    if n > 397
        throw(DomainError(n, "n can't be above 397 or this things barfs") )
        # actually I worry about accuracy before then
    end
    
    d = zeros(n+1) # note indexes offset, so d[1] = d_0
    u = zeros(n+1)
    
    d[1] = 1
    u[1] = 1
    for k=1:n
        u[k+1] = u[k] * 4 * (n+k-1) * (n-k+1)  / ( (2*k) * (2*k-1) )
        d[k+1] = d[k] + u[k+1]
    end
   
    m = 1 .- d ./ d[n+1]

    return m, d, u
end





# Riemann-Siegel algorithm for the critical line for large |t|
#   http://numbers.computation.free.fr/Constants/Miscellaneous/zetaevaluations.html
# but then need to differentiate this :( 



# Also an Euler-Maclaurin formula
#   COMPUTATION OF THE DERIVATIVES OF THE RIEMANN ZETA-FUNCTION IN THE COMPLEX DOMAIN, U.S.S.R. Comput.Maths.Math.Phys.,Vo1.28,No.4,pp.115-124,1988, A.YU. YEREMIN, I.E. KAPORIN and M.K. KERIMOV
#



# need an alternative near 0.0 because
#    reflection doesn't work there because of poles (that cancel) at (1-s) = 1
#    alternating series isn't as accurat
# COMPUTATION OF THE DERIVATIVES OF THE RIEMANN ZETA-FUNCTION IN THE COMPLEX DOMAIN, U.S.S.R. Comput.Maths.Math.Phys.,Vo1.28,No.4,pp.115-124,1988, A.YU. YEREMIN, I.E. KAPORIN and M.K. KERIMOV (equation 3.2)
#   equation 3.1 into 3.2 to avoid the poles
#        see derivatives_zeta.pdf
#
function zeta_d_near_zero( s ; accuracy=1.0e-14)
    # if abs(s) > 1.0e-7
    #     throw(DomainError(s, " assumes abs(s) << 1"))
    # end

    # Taylor series for Cot
    #       https://proofwiki.org/wiki/Power_Series_Expansion_for_Cotangent_Function
    #    cot(x) = 1/x - x/3 + o(x)    # in general cot(x) = sum_{n=0}^\infty B_{2n}/(2n!)  *(2i)^2n * x^(2n-1)
    #    (π/2) * cot(π*s/2) = 1/s - π^2*s/12 - π^4*s^3/720 + ... 
    #
    # Using the Laurent expansion of zeta and its derivative around 1.0 (ie s~0, so 1-s~1)
    #    -ζ'(1-s)/ζ(1-s)
    #            (1/s^2)/(1/s + γ + γ_1 s) 
    #          ~ -1/s - γ - s*(γ^2 + 2 γ_1) - s^2 * ( 3 γ γ_1 + (3/2)*γ_2 + γ^3  )    + O(s^3)  # wolfram alpha
    #                          where γ_i = stieltjes(i)
    # 
    # Taylor series for polygamma:    https://en.wikipedia.org/wiki/Polygamma_function#Taylor_series
    #     ψ(z+1) = -γ + \sum_{k=1}^{\infty} (-1)^{k+1} ζ(k+1) z^k,  for |z|<1
    #    -ψ(1-s) ~ +γ + ζ(2)*s + ζ(3)*s^2
    #               errors here are approx s^3, so works for us for s~ 1.0e-5
    #  

    # note that the γ and the 1/s terms cancel, leaving remainders that behave well around s=0
    # but
    #    cot approximation doesn't seem to be working well
    #    I can't easilty test the zeta approx, but can add terms easily enough
    #  
    
    term1 = log(2π)                  # 
    term2 = zeta(2)*s + zeta(3)*s^2  # the remainder from the polygamma function
    term3 = -π^2* s/12               # the remainder from the cot
    term4 = -s*( γ^2 + 2*stieltjes(1) )  # the remainder from the zetas

    total = zeta(s) * ( term1 + term2 + term3 + term4 )
   
    return total
end


# calculation near the trivial zeros at -2, -4, -6, ...
#   but we don't need to use this now that we have a better reflection formula
#  
function zeta_d_near_zeros( s ; accuracy=1.0e-14)
    m = round( s / 2 )

    if n < 0
        tau = s - 2*m
    else
        error("haven't done these cases yet")
    end

    C = log(2π) - polygamma(0, 1-s) - zeta_derivative( 1 - s; accuracy=accuracy )[1] / zeta(1-s)

    term1 = 2.0 ^(s-1) * pi^s * cos(pi*s/2) * gamma(1-s) * zeta(1-s)
    term2 = 2.0 ^s  * pi^(s-1) * sin(pi*s/2) * gamma(1-s) * zeta(1-s) * C
    
    return term1 + term2, -1000 + n
end


# function to create test values
function test_d_near_zeros( m::Integer )
    # see https://en.wikipedia.org/wiki/Particular_values_of_the_Riemann_zeta_function#Derivatives
    #  "The Riemann Zeta-Function and Its Derivatives", Bejoy K. Choudhury, Choudhury-RiemannZetaFunctionDerivatives-1995.pdf
    #    mainly from OEIS and its sources
    if m==0
        return -log(2π)/2
    elseif m==-2
        return -zeta(3) / (4π^2)
    elseif m==-4
        return 3 * zeta(5) / (4π^4)
    elseif m==-6
        return -45 * zeta(7) / (8π^6)
    elseif m==-8
        return 315 * zeta(9) / (4π^8)
    elseif iseven(m) && m<0
        n = Int64(-m/2)
        total = (-1)^n *  zeta(2*n+1) / 2
        for i=1:2*n
            total *= i/(2π)
        end
        return total
    elseif isodd(m) && m<0 
        # https://www.researchgate.net/profile/Mehmet-Kirdar/publication/362093188_The_Derivative_of_the_Riemann_Zeta_Function_and_Its_Values_at_Integers/links/62d684f47437d7248955a564/The-Derivative-of-the-Riemann-Zeta-Function-and-Its-Values-at-Integers.pdf
        # https://grokipedia.com/page/Particular_values_of_the_Riemann_zeta_function
        n =  Int64(-(m - 1)/2)
        if 2*m > 35
            throw(DomainError(m, "can't cope with very large odd m"))
        end
        term1 = ( polygamma(0, 2*n) - log(2π) ) * bernoulli(2*n) / (2*n)
        term2 = (-1)^(n+1) * factorial(2*n) * zeta_derivative(2*n)[1] / (n*(2π)^(2*n))
        # should change the above into a loop as for even cases
        total = term1 + term2
        return total
    elseif iseven(m) && m > 0
        throw(DomainError(m, "m must be negative"))
        # https://grokipedia.com/page/Particular_values_of_the_Riemann_zeta_function
        #   but this and the odd negative case are circular :( 
        n = Int64(m/2)
        factor = (-1)^(n+1) * (2π)^m / (2 * factorial(m)) 
        term1 = m*test_d_near_zeros( 1 - m )
        term2 = ( polygamma(0,m) - log(2π) ) * bernoulli(m)
        return factor*( term1 - term2 )
    else
        # no closed form for the positive, odd cases
        throw(DomainError(m, "m must be negative"))
    end

    

    # OEIS
    #     zeta'(-n) = (BernoulliB(n+1)*HarmonicNumber(n))/(n+1) - log(A(n))
    #             where
    #                A(n) is the n-th Bendersky constant, that is the n-th generalized Glaisher constant.
    #                
end
