"""
    bernoulli(n, x)

 Calculates Bernoulli polynomials from the Hurwitz-zeta function using

```ζ(-n,x) = -B_{n+1}(x)/(n+1), for Re(x)>0,```

 which is faster than direct calculation of the polynomial except for n<=4, but for 
 small n (n<=5) we use the exact polynomials. For negative x, we use the recursive
 formula to push it into a postive range. 

 e.g., see
 
 + https://en.wikipedia.org/wiki/Bernoulli_polynomials
 + http://dlmf.nist.gov/24

## Arguments
* `n::Integer`: the index into the series, n=0,1,2,3,...
* `x::Real`: the point at which to calculate the polynomial

## Examples
```jldoctest
julia> bernoulli(6, 1.2)
0.008833523809524069
```
"""
function bernoulli(n::Int, x::Real)
    if n<0
        throw(DomainError(n))
    end
    if n == 0
        return 1 # zeta formula doesn't hold for n=0, so return explicit value
    elseif n == 1 # get some easy cases out of the way quickly
        return x-0.5 
    elseif n == 2
        return x^2 - x + 1.0/6.0
    elseif n == 3
        return x^3 - 1.5*x^2 + 0.5*x
    elseif n == 4
        return x^4 - 2.0*x^3 +     x^2 - 1/30.0
    elseif n == 5
         return x^5 - 2.5*x^4 +(5.0/3.0)*x^3 - x/6.0
    end
    # return -n*SpecialFunctions.zeta(1-n, x)
    # see https://carma.newcastle.edu.au/resources/jon/Preprints/Papers/Published-InPress/Oscillatory%20(Tapas%20II)/Papers/coffey-zeta.pdf, p.341
    # however, I am not really happy with the results I am getting out of this
    # so exploit symmetries in B_n(x) to compute recursively for x<=0
    # if x >= 0
    #      return -n*SpecialFunctions.zeta(1-n, x)
    # else
    #     return bernoulli(n, x+1) - n*x^(n-1)
    # end

    # total = 0.0
    # for k = 0:n-1
    #     total -= binomial(n+1, k) *  bernoulli(k,x)
    # end
    # total = total / (n+1)
    # total += x^n
    # return total

    # alternative, a little better, but relies on bernoulli numbers, which I only have to 35
    # and still has some non-negligable errors, so I still would like a better piece in here
    if n <= 35
        total = 0.0
        for k = 0:n
            total += binomial(n, k) *  bernoulli(n-k) * x^k
        end
        return total
    else
        # this is faster but less accurate, and runs into some weird problems, e.g.,
        # SpecialFunctions.zeta(-6, -1.0) = NaN (probably should be 2)
        return -n*SpecialFunctions.zeta(1-n, x)
    end
end
