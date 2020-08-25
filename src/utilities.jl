import Base: parse

"""
 Parse complex numbers.

 This code doesn't deal with all possible forms of complex numbers, just those output by Mathematica

## Arguments
* `::Type{Complex{T}}`: the type to parse to, where T is a Real number type
* `s::AbstractString`: the string to parse

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> parse( Complex{Float64}, "1.2 - 3.1*I")
1.2 + 3.1im
```
"""
function parse(::Type{Complex{T}}, s::AbstractString) where {T<:Real}
    s = replace(s, "*^" => "e") # Mathematica seems to export in forms like "3.061616997868383*^-18"
    # real_r = r"([+-]?\d*(\.\d*)?(e[+-]?\d+)?)([\s+-]|$)"
    real_r = r"([+-]?\d*(\.\d*)?(e[+-]?\d+)?)([^*iI])"
    imag_r = r"([+-]?\d*(\.\d*)?(e[+-]?\d+)?)\s*\*\s*I"
    
    z = zero(Complex{T})
    
    m = match( imag_r, s )
    if m != nothing
        z += im*parse(T, m[1])
    end
    
    m = match( real_r, s )
    if m != nothing
        z += parse(T, m[1])
    end
    return z
end

"""
 Parse rational numbers.

    From https://github.com/JuliaLang/julia/issues/18328

## Arguments
* `::Type{Rational{T}}`: the type to parse to, where T is an Integer number type
* `s::AbstractString`: the string to parse

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> parse( Rational{Int64}, "1 / 2")
1//2
```
"""
function parse(::Type{Rational{T}}, x::AbstractString) where {T<:Integer}
    list = split(x, '/', keepempty=false)
    if length(list) == 1
        return parse(T, list[1]) // 1
    else
        @assert length(list) == 2
        return parse(T, list[1]) // parse(T, list[2])
    end
end
