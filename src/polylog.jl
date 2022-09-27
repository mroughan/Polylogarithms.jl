
const default_accuracy = 1.0e-12
const default_max_iterations = 1000
const near_int_threshold = 1.0e-6
const series_transition_threshold = 0.25

struct Diagnostics # this immutable type has no fields, so constructing it is essentially free
end

"""
    polylog(s, z)

Calculates the Polylogarithm function ``{Li}_s(z)`` defined by
    
``{Li}_s = \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s}``

Uses double precision complex numbers (not arbitrary precision).
It's goal is an relative error bound 10^{-12}.
 
## Input Arguments
* ``s`` `::Complex`: the 'fractional' parameter
* ``z`` `::Complex`: the point at which to calculate it

It should also accept input arguments as `Real` or `Rational` or `Integer` but these aren't completely tested.

There are additional keywords, but these are currently intended for testing not use.

## Output Arguments
* ``Li_s(z)``: The result

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> polylog(0.35, 0.2)
0.23803890574407033
```
"""
function polylog(s::Number, z::Number;
                 level=1, # keep track of recursion
                 accuracy::Float64=default_accuracy,
                 min_iterations::Integer=0,
                 max_iterations::Integer=default_max_iterations)
    polylog(s, z, Diagnostics();
            level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)[1] # just output the result
end

# this is the main version, but outputs diagnosstics, which I guess most people won't want
"""
    polylog(s, z, Diagnostics())

Calculates the Polylogarithm function ``{Li}_s(z)`` defined by

``{Li}_s = \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s}``

Uses double precision complex numbers (not arbitrary precision).
It's goal is an relative error bound 10^{-12}.
 
## Input Arguments
* ``s::`` `Complex`: the 'fractional' parameter
* ``z`` `::Complex`: the point at which to calculate it
* `::Diagnostics`: use this to indicate that the output should include extra information

It should also accept input arguments as `Real` or `Rational` or `Integer` but these aren't completely tested.

There are additional keywords, but these are currently intended for testing not use.

## Output Arguments
* ``Li_s(z)``: The result
* ``n``:       The number of elements used in each series
* `series`:    The series used to compute results (4 = reciprocal)

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> polylog(0.35, 0.2, Diagnostics() )
(0.23803890574407033, 17, 1)
```
"""
function polylog(s::Number, z::Number, ::Diagnostics;
                 level=1, # keep track of recursion
                 accuracy::Float64=default_accuracy,
                 min_iterations::Integer=0,
                 max_iterations::Integer=default_max_iterations)
    tau_threshold = 1.0e-3
    μ = log(convert(Complex{Float64}, z)) # input z could be an integer or anything
    t = abs(μ/twoπ)
    T = 0.512 # the duplication formula seems to work if T=0.5 in the small associated wedge, but wyh risk it?
    if abs(μ) < 1.0e-14
        # deal with this case separately or it causes pain in Series 2 and 3, which should be resolvable, but its just easier here
        # there are lots of other special cases to add here eventually
        if real(s) > 1
            return SpecialFunctions.zeta(s), 0, 0
        else
            return typeof(z)(Inf), 0 , 0
        end
    elseif abs(z) <= 0.5 && abs(z) < t
        return polylog_series_1(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif t <= T && ( abs(round(real(s))-s) > tau_threshold || real(s)<= 0 )
        return polylog_series_2(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif t <= T &&  abs(s) <= tau_threshold # deal with small, positive values of s
        return polylog_series_2(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif t <= T
        return polylog_series_3(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
#    elseif t <= 2.0
    else
        # println("  main level $level, z=$z,  abs(μ)/2π = ", abs(log(z))/twoπ  )
        return polylog_duplication(s, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # elseif abs(z) > 1
    #     return polylog_reciprocal(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # else
    #     throw(DomainError(z, "bad z value $z with |z|=$(abs(z)) and |log(z)/(twoπ)|=$(abs(μ/twoπ))"))
    end
    # we could have a lot more special cases here, particularly for integer input
    # to make the code faster for these cases, but at the moment may use such values for testing
end

# old version without square root
    # if abs(z) <= 0.5 && abs(z) < abs(μ/twoπ)
    #     return polylog_series_1(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # elseif abs(μ/twoπ) < series2_threshold && ( abs(round(real(s))-s) > tau_threshold || real(s)<= 0 )
    #     return polylog_series_2(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # elseif abs(μ/twoπ) < series2_threshold
    #     return polylog_series_3(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # elseif abs(z) > 1
    #     return polylog_reciprocal(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # else
    #     throw(DomainError(z, "bad z value $z with |z|=$(abs(z)) and |log(z)/(twoπ)|=$(abs(μ/twoπ))"))
    # end

# function polylog(s::Number, z::Number, accuracy::Real=default_accuracy)
#     if z ≈ 1.0
#         if real(s) > 1
#             return zeta(s)
#         else
#             return Inf
#         end
#     elseif z ≈ -1.0
#         return -eta(s)
#     elseif s ≈ 0.0
#         return z ./ (1-z)
#     elseif abs(z) <= T
#         ifconfig
#         return polylog_direct(s, z, accuracy)
#     elseif abs(z) >= 1/T && isinteger(s) && real(s) < 0
#         # use reciprocal formula to calculate in terms of Li_n(1/z)
#         # but note for negative integer s, it collapses to something small
#         return -(-1.0)^s .*polylog_direct(s, 1/z, accuracy)
#     elseif  abs(z) >= 1/T
#         # use reciprocal formula to calculate in terms of Li_s(1/z)
#         twopi = 2π
#         z = convert(Complex{Float64}, z)
#         G = (twopi*im)^s * zeta( 1-s, 0.5 + log(-z)/(twopi*im) ) /  gamma(s)
#         F = complex(-1.0)^s * polylog_direct(s, 1/z, accuracy)

#         A = twopi*im*log(z)^(s-1)/(gamma(s))
#         if ( isreal(z) && real(z)>=1 )
#             Θ = 1
#         else
#             Θ = 0
#         end
#         # println("G = $G, F=$F, Θ=$Θ, A=$A")
#         return ( G - F - Θ*A )
#     else 
#         # power series around mu=0, for z = e^mu
#         polylog_series_mu(s, z, accuracy)
#     end
# end

    
####################################################################
#### these are component functions and aren't exported at this point
#### note that for consistency they all have keywords arguments like "accuracy" but
#### these aren't intended for general use, just for testing (at the moment)

# # calculate using the relationship to the Hurwitz zeta function
# function polylog_zeta(s::Number, z::Number, accuracy=default_accuracy)
#     # compute using the Hurwitz-zeta function identity
#     #   N.B. this doesn't seem to work as well as you might think
#     x = im * (log(convert(Complex{Float64}, -z)) / twoπ)
#     ss = 1-s
#     ip = im^ss
#     return ( SpecialFunctions.gamma(ss)/twoπ^(ss) ) * (ip * SpecialFunctions.zeta(ss, 0.5+x) + conj(ip) * SpecialFunctions.zeta(ss, 0.5-x))
# end

# calculate using the duplication formula
function polylog_duplication(s::Number, z::Number;
                             level=0, # keep track of recursion
                             accuracy::Float64=default_accuracy,
                             min_iterations::Integer=0,
                             max_iterations::Integer=default_max_iterations)
    z = convert(Complex{Float64}, z)
    f = min( 0.5, 2.0^(1-real(s)) )
    # println("  dup level $level, z=$z,  abs(μ)/2π = ", abs(log(z))/twoπ  )
    (Li1, k1, series1) = polylog(s,  sqrt(z), Diagnostics();
                                 level=level+1, accuracy=f*accuracy,
                                 min_iterations=min_iterations, max_iterations=max_iterations)
    (Li2, k2, series2) = polylog(s, -sqrt(z), Diagnostics();
                                 level=level+1, accuracy=f*accuracy,
                                 min_iterations=min_iterations, max_iterations=max_iterations)
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into doubles
    end
    return (2^(s-1) * ( Li1 + Li2 ), k1 + k2, 10+series1+series2)
end

# # calculate using the reciprocal formula
# function polylog_reciprocal(s::Number, z::Number;
#                             accuracy::Float64=default_accuracy,
#                             min_iterations::Integer=0,
#                             max_iterations::Integer=default_max_iterations)
#     # z = convert(Complex{Float64}, z)
#     if abs(z) <= 1
#         throw(DomainError(z, "only use this function for |z|>1, and pref |z| > 2"))
#     end
#     # if abs(z) < 2
#     #     warn("Slow convergence for  |z| < 2")
#     # end
#     if abs(s) < 0.1*accuracy
#         return (z/(1-z), 0) # use the identity 
#     elseif real(s) < 0 &&  abs(imag(s)) < 0.1*accuracy && abs( round(s) - s ) < 0.1*accuracy
#         G = 0.0 # pole of the Gamma function
#         A = 0.0
#     else
#         # G = (twoπ*im)^s * SpecialFunctions.zeta( 1-s, 0.5 + log(complex(-z))/(twoπ*im) ) /  SpecialFunctions.gamma(s)
#         # A = twoπ*im*log(z)^(s-1) / SpecialFunctions.gamma(s)
#         tmp = exp( s*log(twoπ*im) - SpecialFunctions.loggamma(complex(s)) )  # (twoπ*im)^s /  SpecialFunctions.gamma(s)
#         G = tmp * SpecialFunctions.zeta( 1-s, 0.5 + log(complex(-z))/(twoπ*im) ) 
#         A = twoπ*im*log(z)^(s-1) / SpecialFunctions.gamma(s)
#     end
#     # accuracy of overall result depends on size of total, which includes these other parts 
#     (Li, k, series) = polylog_series_1(s, 1/z; accuracy=0.1*accuracy,
#                                min_iterations=min_iterations, max_iterations=max_iterations, existing_total=G)
#     F = complex(-1.0)^s * Li 
#     if ( imag(z) == 0 ) &&  ( real(z) >= 1 )
#         Θ = 1.0
#     else 
#         Θ = 0.0
#     end
#     # println("G = $G, F=$F, Θ=$Θ, A=$A")
#     return ( G - F - Θ*A, k, 3+series )
# end

# calculate using direct definition
function polylog_series_1(s::Number, z::Number;
                          accuracy::Float64=default_accuracy,
                          min_iterations::Integer=0,
                          max_iterations::Integer=default_max_iterations,
                          existing_total::Number=0.0)
    # prolly should convert z to a double or complex-double
    if abs(z) > 1 || ( abs(z) ≈ 1  && real(s) <= 2)
        throw(DomainError(z))
    end
    if abs(z) > 1/2
        throw(DomainError(z, "Slow convergence for  |z| > 1/2"))
    end
    total = 0.0
    converged = false
    a = z
    k = 0
    if real(s) < 0
        min_iterations = ceil( real(s) / log(abs(z)) )
    end
    while k<=max_iterations && ~converged
        k = k+1
        total += a
        a *= z * ( k/(k+1.0) )^s
        # println("   total = $total")
        if k > min_iterations && abs(a)/abs(total) < 0.5*accuracy
            converged = true
        end
    end
    return (total, k, 1)
end

# calculate using power series around μ = log(z) = 0
# this should not be used near positive integer values of s, but we allow it here in order to test
function polylog_series_2(s::Number, z::Number;
                          accuracy::Float64=default_accuracy,
                          min_iterations::Integer=0,
                          max_iterations::Integer=default_max_iterations )
    μ = log(convert(Complex{Float64}, z)) # input z could be an integer or anything
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into doubles
    end
    # println("μ = $μ") 
    if abs(μ) > twoπ
        throw(DomainError(z, "we need |log(z)|< 2π for this series"))
    end
    # if real(s) > 0
        # min_iterations = ceil( real(s) ) + 1
    # else
        # min_iterations = ceil( -real(s) ) + 1
        # min_iterations = ceil( real(s) / log(abs( log(z)/twoπ ) ) )
    # end
    oneminuss = 1.0 - s
    total = SpecialFunctions.gamma(oneminuss) * (-μ)^(- oneminuss)
    # total = exp( SpecialFunctions.loggamma(oneminuss) - oneminuss*log(-μ) )
    converged = false
    tmp = 1
    k = 0
    a = Inf
    a_2 = Inf
    # A = abs( 2.0*twoπ^real(s) * exp(abs(imag(π*s)))  )
    # this doesn't work if z=1.0, and hence μ=0, even when that converges, but should already be delt with
    while k<=max_iterations && ~converged
        a_3 = a_2
        a_2 = a
        a = tmp * SpecialFunctions.zeta(s - k)
        total += a
        tmp *= μ/(k+1)
        if k > min_iterations &&
            abs(a)/abs(total) < 0.5*accuracy &&
            abs(a_2)/abs(total) < 0.5*accuracy &&
            abs(a_2) > abs(a)
            # abs( A * (k-real(s))^real(-s) * (μ/twoπ)^k )/abs(total) < accuracy
            # && abs( 2*twoπ^real(s) * (μ/twoπ)^k )/abs(total) < accuracy 
            # the stopping rule should be implemented more efficiently as part of the calculation above
            converged = true
        end
        k = k + 1
    end
    # get correct value along the branch
    if isreal(z) && real(z)>=1
        # total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
        total -= exp( log(twoπ*im) + (s-1)*log(μ) - SpecialFunctions.loggamma(s) )
    end
    return (total, k, 2)
end

function c_closed(n::Integer, j::Integer,  ℒ::Number)
    d2 = SpecialFunctions.digamma(n+1) - ℒ
    if j==0
        return harmonic(n) - ℒ
    elseif j==1
        # Wood:+stieltjes(1) - d2^2/2 + π^2/6 + SpecialFunctions.polygamma(1,n+1)/2 
        return -stieltjes(1) - d2^2/2 - π^2/6 + SpecialFunctions.polygamma(1,n+1)/2 
    elseif j==2
        # Wood:stieltjes(2)   + d2^3/6 + d2*( π^2/6 + SpecialFunctions.polygamma(1,n+1)/2 ) + SpecialFunctions.polygamma(2,n+1)/6
        return stieltjes(2)/2 + d2^3/6 + d2*( π^2/6 - SpecialFunctions.polygamma(1,n+1)/2 ) + SpecialFunctions.polygamma(2,n+1)/6
   end
end

function Q_closed(n::Integer, τ::Number, ℒ::Number; n_terms::Integer=3)
    # τ is the distance from the pole s=n>0, ℒ = log(-μ) = log(-log( z ))
    max_n_terms = 3
    if n_terms < 1 || n_terms > max_n_terms
        throw(DomainError(n_terms))
    end
    return sum( c_closed.(n, 0:n_terms-1,  ℒ) .* τ.^(0:n_terms-1) )
end

function Q(n::Integer, τ::Number, ℒ::Number; n_terms::Integer=5) # Crandall,2012, p.35
    # τ is the distance from the pole s=n>0, ℒ = log(-μ) = log(-log( z ))
    if abs(τ) <= 1.0e-14
        # if really close to the integer, then ignore the extra terms
        return c_closed(n, 0,  ℒ)
    else
        max_n_terms = 7
        if n_terms < 1 || n_terms > max_n_terms
            throw(DomainError(n_terms))
        end
        if n_terms <= 3
            # use the direct method in this case
            return Q_closed(n, τ, ℒ; n_terms=n_terms)
        end
        return sum( c_crandall.(n, 0:n_terms-1,  ℒ) .* τ.^(0:n_terms-1) )
    end
end

function c_crandall(k::Integer, j::Integer,  ℒ) # Crandall,2012, p.35
    return (-1)^j*stieltjes(j)/SpecialFunctions.gamma(j+1) - b_crandall(k,j+1,ℒ)
end

function b_crandall(k::Integer, j::Integer,  ℒ) # Crandall,2012, p.36
    total = 0
    for q=0:j
        for t=0:j-q
            p = j-q-t
            a1 = ℒ^p / SpecialFunctions.gamma(p+1)
            a2 = (-1)^t * f_crandall(k,q) # Bailey and Borwein, 2015 correct Crandall (t+q - > t)
            a3 = g_crandall(t)/SpecialFunctions.gamma(t+1)
            total += a1 * a2 * a3
        end
    end
    return total
end

const gamma_t = [1.0,-0.5772156649015315,1.9781119906559432,-5.44487445648531,23.561474084025583,-117.83940826837748,715.0673625273184,-5019.848872629852,40243.62157333573,-362526.2891146549,3.627042412756892e6,-3.990708415143132e7,4.7894329176518273e8,-6.226641351546061e9,8.717563381070836e10]
function g_crandall(t::Integer) # Crandall,2012, p.17
    # t derivate of Gamma function at 1
    # see "gamma_derivatives.jl" for derivations of these numbers
    if t<0
        throw(DomainError(t))
    elseif t>14
        throw(DomainError(t, "only calculate the 1st 14 derivatives"))
    else
        return gamma_t[t+1]
    end
end

function f_crandall(k::Integer, q::Integer) # Crandall,2012, p.36
    # Crandall omits the k=0 case. but this is given in Bailey and Borwein and matches other text
    if k == 0 && q == 0
        return 1
    elseif k == 0
        return 0
    elseif q == 0
        return 1
    # elseif q == 1
    #     return -harmonic(k)
    # elseif q == 2
    #     return (harmonic(k)^2 + harmonic(k,2))/2
    else
        h = 0:q
        return sum( (-1).^h .* f_crandall.(k-1, q .- h)  ./ k.^h )
    end
end

# For the special case that s is near a positive integer n>0
# Calculate in a power series around z=1, and s=n    
function polylog_series_3(s::Number, z::Number;
                          accuracy::Float64=default_accuracy,
                          min_iterations::Integer=0,
                          max_iterations::Integer=default_max_iterations,
                          n_terms::Integer=5 )
    μ = log(convert(Complex{Float64}, z))
    if abs(μ) > twoπ
        throw(DomainError(z, "does not converge for abs(ln(z)) > twoπ"))
    end
    if real(s)<=0.5
        throw(DomainError(s, "for this function s should be near a positive integer"))
    end
    # this series assumes s is near a positive integer
    n = Int(round(real(s)))
    τ = s - n
    # if real(s) > 0
    #     min_iterations = ceil( real(s) )
    # end
    ℒ = log(complex(-μ))  # '\u2112'
    # total = μ^(n-1)*Q(n-1, τ, ℒ; n_terms=n_terms)/SpecialFunctions.gamma(n)
    total = exp( (n-1)*log(μ) + log(Q(n-1, τ, ℒ; n_terms=n_terms))  - SpecialFunctions.loggamma(n) )
    converged = false
    a = Inf
    a_2 = Inf
    tmp = 1
    k = 0
    while k<=max_iterations && ~converged
        if n - k != 1
            a_2 = a
            a = tmp * SpecialFunctions.zeta(s - k)
            total += a
        end
        tmp *= μ/(k+1)
        if k > min_iterations &&
            abs(a)/abs(total) < 0.5*accuracy &&
            abs(a_2)/abs(total) < 0.5*accuracy &&
            abs(a_2) > abs(a)
            # abs( (μ/twoπ)^k )/abs(total) < 0.05*accuracy
            converged = true
        end
        k = k + 1
    end
    # get correct value along the branch
    if isreal(z) && real(z)>=1 
        total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
    end
    return (total, k, 3)
end

# # asymptotic expansion
# function polylog_series_4(s::Number, z::Number;
#                           accuracy::Float64=default_accuracy,
#                           min_iterations::Integer=0,
#                           max_iterations::Integer=default_max_iterations,
#                           existing_total::Number=0.0)
#     x = log(-convert(Complex{Float64}, z))
#     total = 0.0
#     converged = false
#     k = 0
#     a = Inf
#     while k<=max_iterations && ~converged
#         old_a = a
#         a = (-1)^k * (1.0 - 2.0^(1-2*k)) * twoπ^(2*k) * bernoulli(2*k,0.0) * x^(s-2*k) / ( SpecialFunctions.gamma(2*k+1) * SpecialFunctions.gamma(s+1-2*k) )
#         if abs(a) < abs(old_a)
#             total += a
#         else
#             converged = true
#         end
#         k = k+1
#     end
#     return (total, k, 4)
# end


# ##################################################
    
# # Special case s=integer for use in testing etc. 
#     function polylog(n::Integer, z::Number;
#                       accuracy::Float64=default_accuracy, max_iterations::Integer=default_max_iterations )
#     L = Int(ceil(-log10(accuracy)*log2(10))) # revisit this limit
#     k = 0
#     if n>1
#         μ = log(convert(Complex{Float64}, z))
#         if abs(μ) > twoπ
#             throw(DomainError(z))
#         end
#         total = μ^(n-1)*(harmonic(n-1) - log(-μ))/SpecialFunctions.gamma(n)
#         tmp = 1
#         for m=0:L
#             if n - m != 1
#                 total += tmp * zeta(n - m)
#             end
#             tmp *= μ/(m+1)
#             if abs(tmp)/abs(total) < 1.0e-30
#                 break
#             end
#         end
#         if  isreal(z) && real(z)>=1 
#             total -= 2*pi*im*μ^(s-1)/SpecialFunctions.gamma(n)
#         end
#     elseif n==1
#         total = -log(complex(1-z))
#     elseif n==0
#         total = z / (1-z)
#     elseif n==-1
#         total = z / (1-z)^2
#     elseif n<-1
#         # Crandall's 1.5 for s integer 
#         total = factorial(-n) * (-μ)^(n-1)
#         tmp = 1
#         for k=0:L
#             # total -= μ^k * bernoulli(k-n+1, 0.0) / ( gamma(k+1)*(k-n+1) )
#             total -= tmp * bernoulli(k-n+1, 0.0) / (k-n+1)
#             tmp *= μ/(k+1)
#             if abs(tmp)/abs(total) < 1.0e-30
#                     break
#             end
#         end
#     else
#         error("Should not get this case")
#     end

#     if isreal(z) && real(z)>=1 
#         total -= 2*pi*im*μ^(s-1)/gamma(s)
#     end
#     return (total,k)
# end
