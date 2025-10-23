import Base: parse

"""
 Parse complex numbers.

 This code doesn't deal with all possible forms of complex numbers, just those output by Mathematica

    But, latest version of Mathematica seems to have gone over to a more standard form

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
    real_r = r"([+-]?\d*(\.\d*)?(e[+-]?\d+)?)(?![*iI])"
    # using the negative lookahead to check that it doesn't get followed by a '*', 'i' or an 'I'
    imag_r = r"([+-]?\s*\d*(\.\d*)?(e[+-]?\d+)?)\s*(\*\s*)?[iI]"
    
    z = zero(Complex{T})
    
    m = match( imag_r, s )
    if m != nothing
        tmp = replace( m[1], r"\s*" => "")
        z += im*parse(T, tmp)
    end
    
    m = match( real_r, s )
    if m != nothing
        z += parse(T, m[1])
    end
    return z
end




# drop the next bit as it now seems to be included in Julia's Base
#    see https://github.com/JuliaLang/julia/pull/44550
#    From Jun 17, 2022
#    presuming from version at least v1.10 (the LTS version after this date)
# 
# """
#  Parse rational numbers.
#     From https://github.com/JuliaLang/julia/issues/18328

# ## Arguments
# * `::Type{Rational{T}}`: the type to parse to, where T is an Integer number type
# * `s::AbstractString`: the string to parse

# ## Examples
# ```jldoctest; setup = :(using Polylogarithms)
# julia> parse( Rational{Int64}, "1 / 2")
# 1//2
# ```
# """
# function parse(::Type{Rational{T}}, x::AbstractString) where {T<:Integer}
#     list = split(x, '/', keepempty=false)
#     if length(list) == 1
#         return parse(T, list[1]) // 1
#     else
#         @assert length(list) == 2
#         return parse(T, list[1]) // parse(T, list[2])
#     end
# end
