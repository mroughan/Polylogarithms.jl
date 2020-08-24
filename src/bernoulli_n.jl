"""
    bernoulli(n)

 Calculates the first 34 Bernoulli numbers B_n  (of the first-kind or NIST type) 
 e.g., see

 + http://mathworld.wolfram.com/BernoulliNumber.html
 + https://en.wikipedia.org/wiki/Bernoulli_number
 + http://dlmf.nist.gov/24

 N.B. Bernoulli numbers of second kind only seem to differ in that B_1 = + 1/2 (instead of -1/2)

## Arguments
* `n::Integer`: the index into the series, n=0,1,2,3,...,34 (for larger n use ``bernoulli(n,0.0)`` )

 We only provide the 1st 34 (even) values as beyond this, we can't
 return Int64 rationals, so best to compute the real approximation
 using ``bernoulli(n,0.0). Odd values for n>1 are all zero.``

## Examples
```jldoctest
julia> bernoulli(6)
1 // 42
```
"""
function bernoulli(n::Integer)
    # this just does a lookup -- seemed like it would be easier to code and faster
    # for the size of numbers I am working with
    if n<0
        throw(DomainError(n))
    elseif n > 1 && isodd(n)
        return 0 // 1
    elseif n>34
        throw(DomainError(n, "If n > 34, then the numerator needs Int128 at least, and worse, so this code is not the code you want. Try using bernoulli(n, 0.0) to get a floating point approximation to the result."))
    end

    # Denominator of Bernoulli number B_n
    #   http://oeis.org/A027642
    D = [2, 6, 1, 30, 1, 42, 1, 30, 1, 66, 1, 2730, 1, 6, 1, 510, 1, 798, 1, 330, 1, 138, 1, 2730, 1, 6, 1, 870, 1, 14322, 1, 510, 1, 6, 1, 1919190, 1, 6, 1, 13530, 1, 1806, 1, 690, 1, 282, 1, 46410, 1, 66, 1, 1590, 1, 798, 1, 870, 1, 354, 1, 56786730]

    # Numerator of Bernoulli number B_n (storing 62 of these because they are easy)
    #   http://oeis.org/A027641
    N = [-1, 1, 0, -1, 0, 1, 0, -1, 0, 5, 0, -691, 0, 7, 0, -3617, 0, 43867, 0, -174611, 0, 854513, 0, -236364091, 0, 8553103, 0, -23749461029, 0, 8615841276005, 0, -7709321041217, 0, 2577687858367, 1]
    
    if n==0
        return 1 
    else
        return N[n] // D[n]
    end
end
