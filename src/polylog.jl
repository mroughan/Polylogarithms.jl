
const default_accuracy = 1.0e-12
const default_max_iterations = 1000
const near_int_threshold = 1.0e-6
const series_transition_threshold = 0.25

"""
polylog(s, z)

Calculates the Polylogarithm function ```Li_s(z)``` defined by

    ``` L_s = \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s}```

For ideas going into this see
    
+ Crandall, 'Note on fast polylogarithm computation', 2006, 
  which focusses on the case where s=n (integer and real)    
  http://www.wolfgang-ehrhardt.de/Polylog.pdf

+ Vepstas, 'AN EFFICIENT ALGORITHM FOR ACCELERATING THE CONVERGENCE
  OF OSCILLATORY SERIES, USEFUL FOR COMPUTING THE
  POLYLOGARITHM AND HURWITZ ZETA FUNCTIONS', 2007
  which treats the general case, but presumes arbitrary precision arithmetic
  https://arxiv.org/abs/math/0702243
    
    + Wood, 'The computation of Polylogarithms', 1992
    which focusses on s=n, integer and real, and some typos
    https://www.cs.kent.ac.uk/pubs/1992/110/
    
    + Maximon, 'The dilogarithm function for complex argument', 2003
    which provides useful test cases for s=2.
        
        + Zagier, 'The dilogarithm function in geometry and number theory', 1989 
            similar to Maximon.
            
            Of these the only one that actually specifies a full algorithm is
            Crandall, and he also treats special cases more carefully, so this
            is the one that I have paid most attention to. However, extending it
            for s on the whole complex plane requires some additions, and many
                of these are actually most nicely documented on the wikipedia page
                
                + https://en.wikipedia.org/wiki/Polylogarithm

                With further details at

                + http://mathworld.wolfram.com/Polylogarithm.html
                + http://dlmf.nist.gov/25.12#ii
                + http://mathworld.wolfram.com/Trilogarithm.html
                + http://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/

                The wiki page points out some errors in earlier works, but not all
                parts on the page have references, and not all statements seem to
                come from any of the listed references?

                The code draws heavily on existing functions, in particular the
                Hurwitz-zeta function, which is aliased to zeta(s,q) in Julia.

                    Accuracy has been tested using many of the identities known for Li
                        and relations to known functions as special cases, and by comparison
                        to `polylog(s, z)` in the `mpmath` arbitrary-precision package in Python. 

                        http://mpmath.org/doc/current/functions/zeta.html

                        The latter shows deviations of the order of 

                        + 10^{Im(s) - 20} for Im(s) < 0
                            + 10^{Im(s) - 20} for Im(s) > 0
                                
                                It isn't clear whether we can do better than this with
                                    double-precision arithmetic.

                                    ## Arguments
                                    * `s::Complex`: the 'fractional' coefficient
                                    * `z::Complex`: the point at which to calculate it
                                    * `accuracy::Real=1.0e-18`: nominal accuracy of calculation, but mainly useful for testing

                                        ## Examples
                                        ```jldoctest
    julia> polylog(-1.0, 0.0) 
    (0.0,1)
    ```
                                        """
function polylog(s::Number, z::Number, accuracy::Real=default_accuracy)
    T = 0.25 # threshold at which we change algorithms
    # if !isinteger(s)
    #     return polylog_zeta(s, z)
    # end
    if z ≈ 1.0
        if real(s) > 1
            return zeta(s)
        else
            return Inf
        end
    elseif z ≈ -1.0
        return -eta(s)
    elseif s ≈ 0.0
        return z ./ (1-z)
    elseif abs(z) <= T
        ifconfig
        return polylog_direct(s, z, accuracy)
    elseif abs(z) >= 1/T && isinteger(s) && real(s) < 0
        # use reciprocal formula to calculate in terms of Li_n(1/z)
        # but note for negative integer s, it collapses to something small
        return -(-1.0)^s .*polylog_direct(s, 1/z, accuracy)
    elseif  abs(z) >= 1/T
        # use reciprocal formula to calculate in terms of Li_s(1/z)
        twopi = 2π
        z = convert(Complex{Float64}, z)
        G = (twopi*im)^s * zeta( 1-s, 0.5 + log(-z)/(twopi*im) ) /  gamma(s)
        F = complex(-1.0)^s * polylog_direct(s, 1/z, accuracy)

        A = twopi*im*log(z)^(s-1)/(gamma(s))
        if ( isreal(z) && real(z)>=1 )
            Θ = 1
        else
            Θ = 0
        end
        # println("G = $G, F=$F, Θ=$Θ, A=$A")
        return ( G - F - Θ*A )
    else 
        # power series around mu=0, for z = e^mu
        polylog_series_mu(s, z, accuracy)
    end
end

    
####################################################################
#### these are component functions and aren't exported at this point
#### note that for consistency they all have keywords arguments like "accuracy" but
#### note all make use of these

# calculate using the relationship to the Hurwitz zeta function
function polylog_zeta(s::Number, z::Number, accuracy=default_accuracy)
    # compute using the Hurwitz-zeta function identity
    #   N.B. this doesn't seem to work as well as you might think
    x = im * (log(convert(Complex{Float64}, -z)) / twoπ)
    ss = 1-s
    ip = im^ss
    return ( gamma(ss)/twoπ^(ss) ) * (ip * SpecialFunctions.zeta(ss, 0.5+x) + conj(ip) * SpecialFunctions.zeta(ss, 0.5-x))
end

# calculate using the reciprocal formula
function polylog_reciprocal(s::Number, z::Number;
                            accuracy::Float64=default_accuracy, max_iterations::Int64=default_max_iterations)
    z = convert(Complex{Float64}, z)
    G = (twoπ*im)^s * SpecialFunctions.zeta( 1-s, 0.5 + log(-z)/(twoπ*im) ) /  gamma(s)
    F = complex(-1.0)^s * polylog_series_1(s, 1/z; accuracy=accuracy, max_iterations=max_iterations)
    A = twoπ*im*log(z)^(s-1)/(SpecialFunctions.gamma(s))
    if ( isreal(z) && real(z)>=1 )
        Θ = 1
    else
        Θ = 0
    end
    # println("G = $G, F=$F, Θ=$Θ, A=$A")
    return ( G - F - Θ*A )
end

# calculate using direct definition
function polylog_series_1(s::Number, z::Number;
                          accuracy::Float64=default_accuracy, max_iterations::Int64=default_max_iterations)
    if abs(z) > 1 || ( abs(z) ≈ 1  && real(s) <= 2)
        throw(DomainError(z))
    end
    # if abs(z) > 1/2
    #     warn("Slow convergence for  |z| > 1/2")
    # end
    total = 0.0
    # L = ceil(-log10(accuracy)*log2(10)) # summation limit from Crandall,
    # which is conservative, but based on real s>0
    converged = false
    a = z
    n = 0
    if real(s) < 0
        min_iterations = ceil( real(s) / log(abs(z)) )
    else
        min_iterations = 0
    end
    while n<=max_iterations && ~converged
        n = n+1
        total += a
        a *= z * ( n/(n+1.0) )^s
        # println("   total = $total")
        if n > min_iterations && abs(a)/abs(total) < 0.1*accuracy
            converged = true
        end
    end
    return (total, n)
end

# calculate using power series around μ = log(z) = 0
# this should not be used near positive integer values of s, but we allow it here in order to test
function polylog_series_2(s::Number, z::Number;
                          accuracy::Float64=default_accuracy, max_iterations::Int64=default_max_iterations )
    μ = log(convert(Complex{Float64}, z))
    # println("μ = $μ") 
    if abs(μ) > twoπ
        throw(DomainError(z))
    end
    if real(s) > 0
        min_iterations = ceil( real(s) )
    else
        min_iterations = 0
    end
    total = SpecialFunctions.gamma(1-s) * (-μ)^(s-1)
    converged = false
    tmp = 1
    k = 0
    while k<=max_iterations && ~converged
        a = tmp * SpecialFunctions.zeta(s - k)
        total += a
        tmp *= μ/(k+1)
        if k > min_iterations && abs( (μ/twoπ)^k )/abs(total) < 0.05*accuracy
            # look at the largest value it could have, not including zeta which could drop and then go back up
            converged = true
        end
        k = k + 1
    end
    # get correct value along the branch
    if isreal(z) && real(z)>=1 
        total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
    end
    return (total, k)
end

function c(n::Integer, j::Integer,  ℒ::Number)
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
    return sum( c.(n, 0:n_terms-1,  ℒ) .* τ.^(0:n_terms-1) )
end

# used to be called Q_crandall
function Q(n::Integer, τ::Number, ℒ::Number; n_terms::Integer=5) # Crandall,2012, p.35
    # τ is the distance from the pole s=n>0, ℒ = log(-μ) = log(-log( z ))
    max_n_terms = 5
    if n_terms < 1 || n_terms > max_n_terms
        throw(DomainError(n_terms))
    end
    if n_terms <= 3
        # use the direct method in this case
        return Q_closed(n, τ, ℒ; n_terms=n_terms)
    end
    return sum( c_crandall.(n, 0:n_terms-1,  ℒ) .* τ.^(0:n_terms-1) )
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

function g_crandall(t::Int) # Crandall,2012, p.17
    # t derivate of Gamma function at 1
    #    general form from: http://erikerlandson.github.io/blog/2016/06/15/computing-derivatives-of-the-gamma-function/
    PG21 = -2.4041138063191902 # SpecialFunctions.polygamma(2,1)
    PG41 = -24.8862661234409   # SpecialFunctions.polygamma(4,1)
    if t==0
        return 1
    elseif t==1
        return -γ
    elseif t==2
        # return γ^2 - γ + pi^2/6, Crandall seems to be wrong on this
        return γ^2 + pi^2/6  # from Mathematica
    elseif t==3
        # return -2*γ^3 + 9*γ^2 - (π^2+6)*γ + 3*π^2/2 - 4*SpecialFunctions.zeta(3), Crandall seems to be wrong on this
        return -γ^3  - π^2*γ/2 + PG21 # from Mathematica
    elseif t==4
        return γ^4 + γ^2*π^2 + 3*π^4/20 - 4*γ*PG21 # from Mathematica
    elseif t==5
        return -γ^5 - (20/12)*γ^3*π^2 - (9/12)*γ*π^4 + (10*γ^2 + (20/12)*π^2)*PG21 + PG41 # from Mathematica
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

# For the special case that s is near a postive integer n>0
# Calculate in a power series around z=1, and s=n    
function polylog_series_3(s::Number, z::Number;
                          accuracy::Float64=default_accuracy, max_iterations::Int64=default_max_iterations,
                          n_terms::Integer=5 )
    μ = log(complex(z))
    if abs(μ) > twoπ
        throw(DomainError(z))
    end
    if real(s)<=0.5
        throw(DomainError(s), "for this function s should be near a positive integer")
    end
    # assumes s is near a positive integer
    n = Int(round(real(s)))
    τ = s - n # e = s - n
    if real(s) > 0
        min_iterations = ceil( real(s) )
    else
        min_iterations = 0
    end
    ℒ = log(complex(-μ))  # '\u2112'
    total = μ^(n-1)*Q(n-1, τ, ℒ; n_terms=n_terms)/SpecialFunctions.gamma(n)
    converged = false
    tmp = 1
    k = 0
    while k<=max_iterations && ~converged
        if n - k != 1
            a = tmp * SpecialFunctions.zeta(s - k)
            total += a
        end
        tmp *= μ/(k+1)
        if k > min_iterations && abs( (μ/twoπ)^k )/abs(total) < 0.05*accuracy
            converged = true
        end
        k = k + 1
    end
    # get correct value along the branch
    if isreal(z) && real(z)>=1 
        total -= 2*π*im*μ^(s-1)/SpecialFunctions.gamma(s)
    end
    return (total, k)
end


##################################################
    
# Special case s=integer for use in testing etc. 
    function polylog(n::Integer, z::Number;
                      accuracy::Float64=default_accuracy, max_iterations::Int64=default_max_iterations )
    L = Int(ceil(-log10(accuracy)*log2(10))) # revisit this limit
    k = 0
    if n>1
        μ = log(convert(Complex{Float64}, z))
        if abs(μ) > twoπ
            throw(DomainError(z))
        end
        total = μ^(n-1)*(harmonic(n-1) - log(-μ))/SpecialFunctions.gamma(n)
        tmp = 1
        for m=0:L
            if n - m != 1
                total += tmp * zeta(n - m)
            end
            tmp *= μ/(m+1)
            if abs(tmp)/abs(total) < 1.0e-30
                break
            end
        end
        if  isreal(z) && real(z)>=1 
            total -= 2*pi*im*μ^(s-1)/SpecialFunctions.gamma(n)
        end
    elseif n==1
        total = -log(complex(1-z))
    elseif n==0
        total = z / (1-z)
    elseif n==-1
        total = z / (1-z)^2
    elseif n<-1
        # Crandall's 1.5 for s integer 
        total = factorial(-n) * (-μ)^(n-1)
        tmp = 1
        for k=0:L
            # total -= μ^k * bernoulli(k-n+1, 0.0) / ( gamma(k+1)*(k-n+1) )
            total -= tmp * bernoulli(k-n+1, 0.0) / (k-n+1)
            tmp *= μ/(k+1)
            if abs(tmp)/abs(total) < 1.0e-30
                    break
            end
        end
    else
        error("Should not get this case")
    end

    if isreal(z) && real(z)>=1 
        total -= 2*pi*im*μ^(s-1)/gamma(s)
    end
    return (total,k)
end
