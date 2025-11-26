struct Diagnostics # this immutable type has no fields, so constructing it is essentially free
end

####
#### the following bits and pieces are used for memoization to avoid repeated (expensive) function evaluations
####
struct Memoization # this immutable type has no fields, so constructing it is essentially free
end
import SpecialFunctions.zeta
const CacheZeta  = Dict{Number, ComplexOrReal{Float64} }()
const CacheList = [ CacheZeta ]
# could do the same for gamma (and loggamma, digamma, polygamma, ...) and  

# we need to empty the cache before doing any timing tests
"""
    clearcache()

Note that in this library the zeta function has been replaced by a wrapper that looks up (or stores) the value
into a cache, so that we don't have to repeat zeta function calculations. This function clears the cache so that
you can obtain accurate "first run" performance measurements. It should not be needed in day-to-day calculations unless
it is used a great deal and the cache starts taking up too much memory. 

See also `zeta(s, ::Memoization )`.
"""
function clearcache( ; to_clear = CacheList) 
    # can't empty ExistingSurreals, or it breaks things, and doesn't reduce costs much anyway
    for C in to_clear
        empty!(C)
    end
    return 1
end

# replace zeta function with a version that uses the cache
"""
    zeta(s, ::Memoization )

Note that in this library the zeta function has been replaced by a wrapper that looks up (or stores) the value
into a cache, so that we don't have to repeat zeta function calculations. The traditional zeta function will 
still work, but the code here uses this one, so performance results will be skewed if you don't use the cache correctly. 

See also `clearcache()`.
"""
function zeta( s::Number, ::Memoization  )
    if haskey( CacheZeta, s )
        return CacheZeta[s]
    else
        z = SpecialFunctions.zeta(s)
        CacheZeta[s] = z
        return z
    end
end

####
#### Now start some real code
####
"""
    polylog(s, z)

Calculates the Polylogarithm function ``{Li}_s(z)`` defined by
    
``{Li}_s = \\displaystyle  \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s},``

or by analytic expension to the complex plane. 

It uses double precision complex numbers (not arbitrary precision).
It's goal is an relative error bound 10^{-12}.
 
## Input Arguments
* ``s`` `::Complex`: the 'fractional' parameter
* ``z`` `::Complex`: the point at which to calculate it

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
    
``{Li}_s = \\displaystyle  \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s},``

or by analytic expension to the complex plane. 

It uses double precision complex numbers (not arbitrary precision).
It's goal is an relative error bound 10^{-12}.
 
This version outputs some additional diagnostic information that is useful in debugging, but unlikely to be useful in everyday calculations. 

## Input Arguments
* ``s::`` `Complex`: the 'fractional' parameter
* ``z`` `::Complex`: the point at which to calculate it
* `::Diagnostics`: use this to indicate that the output should include extra information

## Output Arguments
* ``Li_s(z)``: The result
* ``n``:       The number of elements used in each series
* `series`:    The series used to compute results (note this will be a tree when recursion is used
* `max_recursion`:  The maximum depth of recursion used (0 if there is not recursion)

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> polylog(0.35, 0.2, Diagnostics() )
(0.23803890574407033, 17, 1, 0)
```
"""
function polylog(s::Number, z::Number, ::Diagnostics;
                 level=0, # keep track of recursion
                 accuracy::Float64=default_accuracy,
                 min_iterations::Integer=0,
                 max_iterations::Integer=default_max_iterations)
    tau_threshold = 1.0e-3
    μ = log(convert(Complex{Float64}, z)) # input z could be an integer or anything
    t = abs(μ/twoπ)
    T = 0.512 # the duplication formula seems to work if T=0.5 in the small associated wedge, but why risk it?
    # K = 36.84/twoπ # |z|>log(1.0e+16) # conservative bound at which we swap to asymptotic expansion
    s_int = round(real(s)) # nearest real integer to s (presuming it is nearly real)
    if abs(μ) < 1.0e-14
        # Deal with this case separately or it causes pain in Series 2 and 3, which should be resolvable,
        # but its just easier here.
        if real(s) > 1
            return zeta(s, Memoization() ), 0, 0, 0
        else
            return typeof(z)(Inf), 0, 0, 0
        end
        # There are lots of other special cases to add here eventually
    elseif abs(z) <= 0.5 && abs(z) < t
        return polylog_series_1(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif t <= T && ( abs(s_int-s) > tau_threshold || real(s)<= 0 )
        return polylog_series_2(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif t <= T &&  abs(s) <= tau_threshold # deal with small, positive values of s
        return polylog_series_2(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    # elseif t > K 
    #     return  polylog_asympt_series_1(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif t <= T
        return polylog_series_3(s, z; accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    elseif abs(s_int-s) > tau_threshold || real(s)<= 0
        return polylog_root(s, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    else
        # use duplication for now, but should prolly be replaced by a mth-root-m of series 3
        return polylog_duplication(s, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
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
# #   N.B. this doesn't seem to work as well as you might think
# # Julia's "generalised Hurwitz zeta" is the second form (called Zeta in Mathematica), not the one we need
# function polylog_zeta(s::Number, z::Number, accuracy=default_accuracy)
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
#     f = 2.0^(-real(s)) # does bad things
#     f = 1.0 # does bad things
    # println("  dup level $level, z=$z,  abs(μ)/2π = ", abs(log(z))/twoπ  )
    (Li1, k1, series1, max_recursion1) = polylog(s,  sqrt(z), Diagnostics();
                                                 level=level+1, accuracy=f*accuracy,
                                                 min_iterations=min_iterations,
                                                 max_iterations=max_iterations)
    (Li2, k2, series2, max_recursion2) = polylog(s, -sqrt(z), Diagnostics();
                                                 level=level+1, accuracy=f*accuracy,
                                                 min_iterations=min_iterations,
                                                 max_iterations=max_iterations)
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into a (complex) double
    end
    max_recursion = max( level, max_recursion1, max_recursion2 )
    return (2^(s-1) * ( Li1 + Li2 ), k1 + k2, (series1, series2), max_recursion )
end

# calculate using the m-th root formula, naive version
#   this doesn't fix the cancellation problem of the duplication formula
#   so it is only included as a reference
function polylog_root_old(s::Number, z::Number;
                          level=0, # keep track of recursion
                          accuracy::Float64=default_accuracy,
                          min_iterations::Integer=0,
                          max_iterations::Integer=default_max_iterations,
                          α = 1.2)
    # α = 1.2 inconvenient compromise, could result in one step of duplication formula,
    # but if α is too close to 1.0 it inflates m a lot
    z = convert(Complex{Float64}, z)
    μ = log(z)
    m = ceil( real(μ) / (sqrt(α^2-1)*pi) )
    f = min( 0.5, m^(1-real(s)) )
    zm = exp(μ/m)
    L = 0
    k = 0
    series = []
    max_recursion = level
    for j = 0 : 1 : m - 1
        x = zm * exp( (im * 2 * pi * j )/ m )
        (Li1, k1, series1, max_recursion1) = polylog(s,  x, Diagnostics();
                                                     level=level+1, accuracy=f*accuracy,
                                                     min_iterations=min_iterations,
                                                     max_iterations=max_iterations)
        L += Li1
        k += k1
        series = push!( series, series1 )
        max_recursion = max( max_recursion, max_recursion1)
    end
    
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into a (complex) double
    end
    return ( m^(s-1) * L, k, Tuple(series), max_recursion )
end

# calculate using the m-th root formula, summation order reversal in series 2
#   presuming s isn't near integer values ...
#
# also note that at somewhere around 1.0e62, we end up with m=~69, and things start breaking down
#  
function polylog_root(s::Number, z::Number;
                      level=0, # keep track of recursion
                      accuracy::Float64=default_accuracy,
                      min_iterations::Integer=0,
                      max_iterations::Integer=default_max_iterations,
                      α = 1.2)
    # α = 1.2 inconvenient compromise, could result in one step of duplication formula,
    # but if α is too close to 1.0 it inflates m a lot
    
    if typeof(s) <: Real
        s = convert(Float64, s) # convert s into a double
    elseif typeof(s) <: Complex
        s = convert(Complex{Float64}, s) # convert s into a (complex) double
    end
    z = convert(Complex{Float64}, z)
    μ = log(z)
    m = max( 2, Integer(ceil( real(μ) / (sqrt(α^2-1)*pi) )) ) # no point in this series unless m>1
    μ_m = μ/m
    k = 0
    series = []
    max_recursion = level
    ell = Integer(ceil( -m/2 - imag(μ)/twoπ ))
                                      # ell such that
                                      #   -π <=  (imag(μ) + 2*π*ell)/m
                                      #          (imag(μ) + 2*π*(ell+m-1))/m < π

    # first term
    oneminuss = 1.0 - s
    eval_points = [ μ_m  +  twoπ * im * (j/m)  for  j = ell : ell + m - 1 ]
    for i=1:m
        tmp = eval_points[i]
        if imag(tmp) < -π 
            error("  bad choice (1) of ℓ = $ell, m=$m, i=$i, μ=$μ, μ_m=$μ_m, tmp=$tmp")
        end
        if imag(tmp) >= π
            error("  bad choice (2) of ℓ = $ell, m=$m, i=$i, μ=$μ, μ_m=$μ_m, tmp=$tmp")
        end
    end
    total = 0.0
    for i=1:m
        total += (-eval_points[i])^(-oneminuss)
    end
    total *= SpecialFunctions.gamma(oneminuss)

    # series
    converged = false
    tmp = 1
    k = 0
    a = Inf
    a_2 = Inf
    summation_terms = ones(ComplexF64, m)
    summation_terms_old = ones(ComplexF64, m)
    while k<=max_iterations && ~converged
        a_2 = a
        a = sum(summation_terms)
        a *= zeta(s - k, Memoization())
        total += a

        for i=1:m
            tmp = eval_points[i]
            summation_terms_old[i] = summation_terms[i]
            summation_terms[i] *= tmp / (k+1)
        end
        
        if k > min_iterations &&
               abs(a)/abs(total) < 0.5*accuracy &&
               abs(a_2)/abs(total) < 0.5*accuracy &&
               abs(a_2) > abs(a)
            converged = true
        end
        k = k + 1
    end
    L = m^(s-1) * total
    
    # get correct value along the branch
    if isreal(z) && real(z)>=1
        # total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
        L -= exp( log(twoπ*im) + (s-1)*log(μ) - SpecialFunctions.loggamma(s) )
    end
 
    series = 2 + 10*m # use this to signal the value of m used with series 2
    return (L, k, series, max_recursion )
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
    series = 1
    max_recursion = 0
    return (total, k, series, max_recursion)
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
        a = tmp * zeta(s - k, Memoization())
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
    
    series = 2
    max_recursion = 0
    return (total, k, series, max_recursion)
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
            a = tmp * zeta(s - k, Memoization() )
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
    series = 3
    max_recursion = 0
    return (total, k, series, max_recursion)
end

#
# asymptotic expansion
#     
function polylog_asympt_series_1(s::Number, z::Number;
                                 accuracy::Float64=default_accuracy,
                                 min_iterations::Integer=0,
                                 max_iterations::Integer=default_max_iterations,
                                 existing_total::Number=0.0)
    x = log(-convert(Complex{Float64}, z))
    total = 0.0
    converged = false
    k = 0
    a = Inf
    while k<=max_iterations && ~converged
        old_a = a
        a = (-1)^k * (1.0 - 2.0^(1-2*k)) * twoπ^(2*k) * bernoulli( Int128(2*k) ) * x^(s-2*k) /
                  ( SpecialFunctions.gamma(2*k+1) * SpecialFunctions.gamma(s+1-2*k) )
        if k>min_iterations &&
            (abs(a) >= abs(old_a) || abs(a)/abs(total) < 0.5*accuracy || k>max_iterations )
            converged = true
        else
            total += a
            k = k+1
        end
    end
    # get correct value along the branch
    if isreal(z) && real(z)>=1 
        μ = log(convert(Complex{Float64}, z))
        total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
    end 
    series = 4
    max_recursion = 0
    return (total, k, series, max_recursion)
end

 
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
 

##################
### derivatives
##################



"""
    polylog_dz(s, z)

Derivative of the Polylogarithm function ``{Li}_s(z)`` with respect to z.

Note that (see eg [https://en.wikipedia.org/wiki/Polylogarithm](https://en.wikipedia.org/wiki/Polylogarithm))
`` \\frac{d}{dz} Li_s(z) = Li_{s-1}(z)/z ``
   
## Input Arguments
* ``s`` `::Complex`: the 'fractional' parameter
* ``z`` `::Complex`: the point at which to calculate it

## Output Arguments
* ``\\displaystyle \\frac{d}{dz} Li_s(z)``

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> polylog_dz(0.35, 0.2)
1.421095587670745
```
"""
function polylog_dz(s::Number, z::Number;
                    level=1, # keep track of recursion
                    accuracy::Float64=default_accuracy,
                    min_iterations::Integer=0,
                    max_iterations::Integer=default_max_iterations)
    return polylog(s-1, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations) / z
end



"""
    polylog_ds(s, z)

Derivative of the Polylogarithm function ``{Li}_s(z)`` with respect to s.

Note that this is not replicating all the work above (yet), and only doing the
simple series version which is valid only for |z| < 1.
   
## Input Arguments
* ``s`` `::Complex`: the 'fractional' parameter
* ``z`` `::Complex`: the point at which to calculate it

## Output Arguments
* ``\\displaystyle \\frac{d}{ds} Li_s(z)``

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> polylog_ds(0.35, 0.2)
-0.02947228342617501
```
"""
function polylog_ds(s::Number, z::Number;
                    level=1, # keep track of recursion
                    accuracy::Float64=default_accuracy,
                    min_iterations::Integer=0,
                    max_iterations::Integer=default_max_iterations)
    
    if abs(z) >= 1.0
        throw(DomainError(z, "At present, this only works for |z|<1"))
    end
    return polylog_ds_series_1(s, z;
                               accuracy=accuracy,
                               min_iterations=min_iterations,
                               max_iterations=max_iterations)[1]
end

# calculate using direct definition
function polylog_ds_series_1(s::Number, z::Number;
                             accuracy::Float64=default_accuracy,
                             min_iterations::Integer=0,
                             max_iterations::Integer=100000000) # temporarily have made this large because we are pushing towards z->1.0
    total = 0.0
    converged = false
    k = 2
    a = z^k * log(k) / k^s
    if real(s) < 0
        min_iterations = ceil( real(s) / log(abs(z)) )
    end
    while k<=max_iterations && ~converged
        total -= a
        k = k+1
        a *= z * ( (k-1)/(k) )^s * log(k)/log(k-1)
        # println("   total = $total")
        if k > min_iterations && abs(a)/abs(total) < 0.5*accuracy
            converged = true
        end
    end

    # k = 2
    # while k<=max_iterations && ~converged
    #     a = z^k * log(k) / k^s
    #     total -= a
    #     k = k+1
    #     # println("   total = $total")
    #     if k > min_iterations && abs(a)/abs(total) < 0.5*accuracy
    #         converged = true
    #     end
    # end
    
    series = 101 # derivative series 1
    max_recursion = 0
    return (total, k, series, max_recursion)
end
