
default_accuracy = 1.0e-12;

"""
    polylog(s, z)

 Calculates the Polylogarithm function ``Li_s(z)`` defined by

`` L_s = \\sum_{n=1}^{\\infty} \\frac{z^n}{n^s}``

 For ideas going into this see
       
 + Crandall, "Note on fast polylogarithm computation", 2006, 
           which focusses on the case where s=n (integer and real)    
           http://www.wolfgang-ehrhardt.de/Polylog.pdf

 + Vepstas, "AN EFFICIENT ALGORITHM FOR ACCELERATING THE CONVERGENCE
                 OF OSCILLATORY SERIES, USEFUL FOR COMPUTING THE
                 POLYLOGARITHM AND HURWITZ ZETA FUNCTIONS", 2007
           which treats the general case, but presumes arbitrary precision arithmetic
           https://arxiv.org/abs/math/0702243
             
 + Wood, "The computation of Polylogarithms", 1992
           which focusses on s=n, integer and real, and which has formatting issues making
           it hard to read correctly.
           https://www.cs.kent.ac.uk/pubs/1992/110/

 + Maximon, "The dilogarithm function for complex argument", 2003
           which provides useful test cases for s=2

 + Zagier,  "The dilogarithm function in geometry and number theory", 1989 
           similar to Maximon
           
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
    x = im * (log(complex(-z)) / twoπ)
    ss = 1-s
    ip = im^ss
    return ( gamma(ss)/twoπ^(ss) ) * (ip * zeta(ss, 0.5+x) + conj(ip) * zeta(ss, 0.5-x))
end

    
# calculate using direct definition
function polylog_series_1(s::Number, z::Number, accuracy=1.0e-15)
    if abs(z) > 1 || ( abs(z) ≈ 1  && real(s) <= 2)
        throw(DomainError(z))
    end
    if abs(z) > 1/2
        warn("Slow convergence for  |z| > 1/2")
    end
    total = 0.0
    L = ceil(-log10(accuracy)*log2(10)) # summation limit from Crandall, which is conservative, but based on real s
    a = z;
    for n=1:L
        total += a
        a *= z * ( n/(n+1.0) )^s
        # println("   total = $total")
        if abs(a)/abs(total) < 1.0e-30
            break
        end
    end
    return total
end

# calculate using power series around μ = log(z) = 0
function polylog_series_2(s::Number, z::Number, accuracy=default_accuracy)
   μ = log(convert(Complex{Float64}, z))
    # println("μ = $μ") 
    if abs(μ) > 2*pi
        throw(DomainError(z))
    end
    near_int_threshold = 1.0e-4 # nearer than this to the s=in>0, and we need to do something special
    L = Int(ceil(-log10(accuracy)*log2(10))) # revisit this limit
    if isinteger(s)
        n = Int(round(s))
        if n>1
            # Crandall's 1.4 for s integer
            total = μ^(n-1)*(harmonic(n-1) - log(-μ))/gamma(n)
            # println("   μ=$μ, total = $total")
            tmp = 1
            for m=0:L
                if n - m != 1
                    # total += μ^m * zeta(n - m) / gamma(m+1)
                    total += tmp * zeta(n - m)
                end
                # println("   m=$m, total = $total, tmp=$tmp, ctmp=$(μ^m /gamma(m+1))")
                tmp *= μ/(m+1)
                if abs(tmp)/abs(total) < 1.0e-30
                    break
                end
            end
            # println("   μ=$μ, total = $total")
            if  isreal(z) && real(z)>=1 
                total -= 2*pi*im*μ^(s-1)/gamma(n)
            end
            # println("   μ=$μ, total = $total")
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

    elseif real(s) > 1 && abs(s - round(real(s))) < near_int_threshold
        # should have a case in here for when s is close to a real, positive integer (see Wood 9.4, but note it is wrong)
        # basically, we take out the two terms that individually go to inf, but whose sum doesn't
        # and which in the integer case gives μ^(n-1)*(harmonic(n-1) - log(-μ))/gamma(n)
        # and we expand around these
        n = Int(round(real(s)))
        if n==1
            hn = 0
        else
            hn = harmonic(n-1)
        end
        e = s - n
        lmu = log(-μ)
        lmu2 = log(μ)sti
        Q1 = (hn - lmu)
        # wood's version
        # Q2 =  -stieltjes(1) - (digamma(n)-lmu)^2/2 + pi^2/6 + trigamma(n)/2  # digamma(x) = polygamma(0,x), trigamma(x) = polygamma(1,x)
        Q2 = -stieltjes(1) - (digamma(n)-lmu)^2/2 - pi^2/6 + polygamma(1,n)/2 


        # which of these (prolly second)
        # Q3 =  stieltjes(2) + (digamma(n)-lmu)^3/6 + (pi^2/6 + polygamma(1,n)/2)*(digamma(n)-lmu) + polygamma(2,n)/2
        Q3 =  stieltjes(2) + (digamma(n)-lmu)^3/6 + (pi^2/6 + polygamma(1,n)/2)*(digamma(n)-lmu) + polygamma(2,n)/6

        
        Q = Q1 + e*Q2 * e^2*Q3
        # Q = Q1 + e*Q2 
        total = μ^(n-1)*Q/gamma(n)

        # println("   μ=$μ, total = $total")
        tmp = 1
        for k=0:L
            if n - k != 1
                # total += μ^m * zeta(n - k) / gamma(k+1)
                total += tmp * zeta(s - k)
            end
            # println("   m=$m, total = $total, tmp=$tmp, ctmp=$(μ^m /gamma(m+1))")
            tmp *= μ/(k+1)
            if abs(tmp)/abs(total) < 1.0e-30
                break
            end
        end
        # println("   μ=$μ, total = $total") μ = log(convert(Complex{Float64}, z))
        if  isreal(z) && real(z)>=1 
            total -= 2*pi*im*μ^(s-1)/gamma(s)
        end
        # println("   total = $total")
    else
        # equivalent of Crandalls 1.4 for s non-integer, Wood (9.3)
        total = gamma(1-s) * (-μ)^(s-1)
        # println("   μ=$μ, total = $total")
        tmp = 1
        for k=0:L
            # total += μ^k * zeta(s-k)/factorial(Float64(k))
            total += tmp * zeta(s - k)
            # println("      tmp=$(tmp* zeta(s-k)),  total = $total")
            tmp *= μ/(k+1)
            if abs(tmp)/abs(total) < 1.0e-30
                break
            end
        end
        if isreal(z) && real(z)>=1 
            total -= 2*pi*im*μ^(s-1)/gamma(s)
        end
    end
    return total
end

# calculate using power series around μ = log(z) = 0
# for the special case that s=n>0
function polylog_series_3(s::Number, z::Number, accuracy=default_accuracy)
    μ = log(convert(Complex{Float64}, z))
    # println("μ = $μ") 
    if abs(μ) > 2*pi
        throw(DomainError(z))
    end
    near_int_threshold = 1.0e-4 # nearer than this to the s=in>0, and we need to do something special
    L = Int(ceil(-log10(accuracy)*log2(10))) # revisit this limit
    if isinteger(s)
        n = Int(round(s))
        if n>1
            # Crandall's 1.4 for s integer
            total = μ^(n-1)*(harmonic(n-1) - log(-μ))/gamma(n)
            # println("   μ=$μ, total = $total")
            tmp = 1
            for m=0:L
                if n - m != 1
                    # total += μ^m * zeta(n - m) / gamma(m+1)
                    total += tmp * zeta(n - m)
                end
                # println("   m=$m, total = $total, tmp=$tmp, ctmp=$(μ^m /gamma(m+1))")
                tmp *= μ/(m+1)
                if abs(tmp)/abs(total) < 1.0e-30
                    break
                end
            end
            # println("   μ=$μ, total = $total")
            if  isreal(z) && real(z)>=1 
                total -= 2*pi*im*μ^(s-1)/gamma(n)
            end
            # println("   μ=$μ, total = $total")
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

    elseif real(s) > 1 && abs(s - round(real(s))) < near_int_threshold
        # should have a case in here for when s is close to a real, positive integer (see Wood 9.4, but note it is wrong)
        # basically, we take out the two terms that individually go to inf, but whose sum doesn't
        # and which in the integer case gives μ^(n-1)*(harmonic(n-1) - log(-μ))/gamma(n)
        # and we expand around these
        n = Int(round(real(s)))
        if n==1
            hn = 0
        else
            hn = harmonic(n-1)
        end
        e = s - n
        lmu = log(-μ)
        lmu2 = log(μ)sti
        Q1 = (hn - lmu)
        # wood's version
        # Q2 =  -stieltjes(1) - (digamma(n)-lmu)^2/2 + pi^2/6 + trigamma(n)/2  # digamma(x) = polygamma(0,x), trigamma(x) = polygamma(1,x)
        Q2 = -stieltjes(1) - (digamma(n)-lmu)^2/2 - pi^2/6 + polygamma(1,n)/2 


        # which of these (prolly second)
        # Q3 =  stieltjes(2) + (digamma(n)-lmu)^3/6 + (pi^2/6 + polygamma(1,n)/2)*(digamma(n)-lmu) + polygamma(2,n)/2
        Q3 =  stieltjes(2) + (digamma(n)-lmu)^3/6 + (pi^2/6 + polygamma(1,n)/2)*(digamma(n)-lmu) + polygamma(2,n)/6

        
        Q = Q1 + e*Q2 * e^2*Q3
        # Q = Q1 + e*Q2 
        total = μ^(n-1)*Q/gamma(n)

        # println("   μ=$μ, total = $total")
        tmp = 1
        for k=0:L
            if n - k != 1
                # total += μ^m * zeta(n - k) / gamma(k+1)
                total += tmp * zeta(s - k)
            end
            # println("   m=$m, total = $total, tmp=$tmp, ctmp=$(μ^m /gamma(m+1))")
            tmp *= μ/(k+1)
            if abs(tmp)/abs(total) < 1.0e-30
                break
            end
        end
        # println("   μ=$μ, total = $total") μ = log(convert(Complex{Float64}, z))
        if  isreal(z) && real(z)>=1 
            total -= 2*pi*im*μ^(s-1)/gamma(s)
        end
        # println("   total = $total")
    else
        # equivalent of Crandalls 1.4 for s non-integer, Wood (9.3)
        total = gamma(1-s) * (-μ)^(s-1)
        # println("   μ=$μ, total = $total")
        tmp = 1
        for k=0:L
            # total += μ^k * zeta(s-k)/factorial(Float64(k))
            total += tmp * zeta(s - k)
            # println("      tmp=$(tmp* zeta(s-k)),  total = $total")
            tmp *= μ/(k+1)
            if abs(tmp)/abs(total) < 1.0e-30
                break
            end
        end
        if isreal(z) && real(z)>=1 
            total -= 2*pi*im*μ^(s-1)/gamma(s)
        end
    end
    return total
end
