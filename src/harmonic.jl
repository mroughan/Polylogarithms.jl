"""
    harmonic(n::Integer)

 Calculates harmonic numbers,
   e.g., see [http://mathworld.wolfram.com/HarmonicNumber.html](http://mathworld.wolfram.com/HarmonicNumber.html)

## Arguments
* ``n`` `::Integer`: non-negative index of the Harmonic number to calculate

## Examples
```jldoctest
julia> harmonic(2)
1.5
```
"""
function harmonic(n::Integer)
    if n < 0
        throw(DomainError(n))
    elseif n==0
        return 0.0 # by convention (see Crandall, p.22)
    elseif n <= 10
        # perform exact sum for small n
        total = 0.0
        for k=1:n
            total +=  1.0 / k
        end
        return total
    else
        return γ + SpecialFunctions.digamma(n+1) # digamma(m) = ψ(m)
    end
end

"""
    harmonic(x::ComplexOrReal{Float64})

 Calculates harmonic numbers extended to non-integer arguments using the
 digamma form.

## Arguments
* ``x`` `::ComplexOrReal{Float64}`: index of the Harmonic number to calculate

## Examples
```jldoctest
julia> harmonic(2.0)
1.5000000000000016
```
"""
function harmonic(x::ComplexOrReal{Float64})
    return γ + SpecialFunctions.digamma(x+1)
end
#        ``digamma(m,x) = (-1)^{m+1} m! \zeta(m+1,x), m>=1``, and integer

"""
    harmonic(n::Integer,r::Real)

 Calculates generalized harmonic numbers, 
   e.g., see [http://mathworld.wolfram.com/HarmonicNumber.html](http://mathworld.wolfram.com/HarmonicNumber.html)

## Arguments
* ``n`` `::Integer`: non-negative index 1 of the Harmonic number to calculate
* ``r`` `::Real`: index 2 of the Harmonic number to calculate

It should be possible to extend this to complex r, but that requires more testing.

## Examples
```jldoctest
julia> harmonic(2,1.5)
1.3535533905932737
```
"""
function harmonic(n::Integer, r::Real)
    if n < 0
        throw(DomainError(n))
    end
    if n == 0
        return 0.0
    end
    if r==1
        return harmonic(n)
    end
    total = 0.0
    for k=1:n
        total +=  1.0 / k^r
    end
    return total
end


"""
    harmonic(n::Integer,r::Integer)

 Calculates generalized harmonic numbers
   e.g., see [http://mathworld.wolfram.com/HarmonicNumber.html](http://mathworld.wolfram.com/HarmonicNumber.html)
 using a better approach which works when both inputs are integers
 [https://carma.newcastle.edu.au/resources/jon/Preprints/Papers/Published-InPress/Oscillatory%20(Tapas%20II)/Papers/coffey-zeta.pdf](https://carma.newcastle.edu.au/resources/jon/Preprints/Papers/Published-InPress/Oscillatory%20(Tapas%20II)/Papers/coffey-zeta.pdf), p.341
 
## Arguments
* ``n`` `::Integer`: non-negative index 1 of the Harmonic number to calculate
* ``r`` `::Integer`: index 2 of the Harmonic number to calculate

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> harmonic(2,1)
1.5000000000000002
```
"""
function harmonic(n::Integer, r::Integer)
    if r<1
        throw(DomainError(r))
    end
    return (-1)^(r-1) * ( SpecialFunctions.polygamma(r-1,n+1) - SpecialFunctions.polygamma(r-1,1) ) / SpecialFunctions.gamma(r)
end

# # better approach to calculation, but doesn't seem to work
# #  using 
# function harmonic3(x::Real, r::Integer)
#     if (r < 1)
#         throw(DomainError(r))
#     end
#     return SpecialFunctions.zeta(r,1) - SpecialFunctions.zeta(r,x+1)
# end

