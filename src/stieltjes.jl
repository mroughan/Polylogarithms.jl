"""
    stieltjes(n)

 Provides the first 10 Stieltjes (generalized Euler-Mascheroni) constants (see
 Abramowitz and Stegunm, 23.2.5) or [https://en.wikipedia.org/wiki/Stieltjes_constants](https://en.wikipedia.org/wiki/Stieltjes_constants).

 There is a table at "The Generalized Euler-Mascheroni Constants", O.R. Ainsworth and L.W.Howell
 NASA Technical Paper 2264, Jan 1984 
[https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19840007812.pdf](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19840007812.pdf)
 but the OEIS has more accurate values, which will be useful when I get around to 
 bit float versions of the code.

 Note that stieltjes(0) = γ, the Euler–Mascheroni constant, also called just Euler's constant. 
 [https://en.wikipedia.org/wiki/Euler-Mascheroni_constant](https://en.wikipedia.org/wiki/Euler%E2%80%93Mascheroni_constant)

## Arguments
* `n::Integer`: the number of elements to compute.

## Examples
```jldoctest
julia> stieltjes(0)
0.5772156649015329
```
"""
function stieltjes(n::Integer)
    if n<0 || n>10
        throw(DomainError(n))
    elseif n>10
        throw(DomainError(n, "Only the fest 11 Stieltjes numbers are defined so far."))
    end
    # old values from Ainsworth and Howell
    # gamma = [0.57721566490153
    #          -0.07281584548368,
    #          -0.00969036319287,
    #          0.00205383442030,
    #          0.00232537006546,
    #          0.00079332381728,
    #          -0.0002387693455,
    #          -0.0005272895671,
    #          -0.0003521233539,
    #          -0.0000343947747,
    #          ]

    # at the moment these are just Float64, so most of the precision below is lost
    # would be better to make these "Irrational" but that can wait a little as there
    # are worse pieces of the code
    stieltjes = [γ, 
             -0.0728158454836767248605863758749013191377363383, 	# A082633
             -0.0096903631928723184845303860352125293590658061, 	# A086279
             +0.0020538344203033458661600465427533842857158044, 	# A086280
             +0.0023253700654673000574681701775260680009044694, 	# A086281
             +0.0007933238173010627017533348774444448307315394, 	# A086282
             -0.0002387693454301996098724218419080042777837151, 	# A183141
             -0.0005272895670577510460740975054788582819962534, 	# A183167
             -0.0003521233538030395096020521650012087417291805, 	# A183206
             -0.0000343947744180880481779146237982273906207895, 	# A184853
             +0.0002053328149090647946837222892370653029598537, 	# A184854
             ]
    return stieltjes[n+1]
end

# old crude computation
# function gen_euler_calc()
#     # N.B. by default
#     #  log^0(1) = 1
#     for n=0:10
#         m = 100000000 # just choose a big value for the moment
#         total = -log(m)^(n+1) / (n+1)
#         for k=1:m
#             total += log(k)^n / k
#         end
#     end
# end
