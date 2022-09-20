# The Bernoulli numbers +
#    throw in the Euler numbers for free

"""
    bernoulli(n)

 Provides a lookup table for the first 59 Bernoulli numbers ``B_n``  (of the first-kind or NIST type) 
 e.g., see

 + [http://mathworld.wolfram.com/BernoulliNumber.html](http://mathworld.wolfram.com/BernoulliNumber.html)
 + [https://en.wikipedia.org/wiki/Bernoulli_number](https://en.wikipedia.org/wiki/Bernoulli_number)
 + [http://dlmf.nist.gov/24](http://dlmf.nist.gov/24)
 + Abramowitz and Stegun, Handbook of Mathematical Function, ..., Table 23.2

 N.B. Bernoulli numbers of second kind only seem to differ in that ``B_1 = + 1/2`` (instead of -1/2)

 There is an algorithmic version of this, but given that the number which are reasonable sized is small, a lookup tables seems fastest. 

## Arguments
* ``n`` `::Integer`: the index into the series, ``n=0,1,2,3,...,59`` (for larger ``n`` use `bernoulli(n,0.0)` )

 The type of the output is Rational{typeof(n)}, so we can only calculate number such that this would not cause
 round-off, i.e,
     typeof(n)     max_n
        Int32       23
        Int64       35
        Int128      59

 Beyond this, its probably best to compute the real approximation
 using `bernoulli(n,0.0)`. 

 Harvey has an algorithm used to get n=100,000,000 but this seems overkill for what I need. 
   + Harvey, David (2010), "A multimodular algorithm for computing Bernoulli numbers", Math. Comput., 79 (272): 2361–2370, arXiv:0807.1347, doi:10.1090/S0025-5718-2010-02367-1, S2CID 11329343, Zbl 1215.11016
   + Apparently implemented in SageMath (since 3.1)

 But odd values for ``n>1`` are all zero, so they are easy.

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> bernoulli(6)
1//42
```
"""
function bernoulli(n::Integer)
    if n<0 
        throw(DomainError(n))
    elseif n > 1 && isodd(n)
        return 0 // 1
    elseif typeof(n)==Int32 && n > 23
        throw(DomainError(n, "If n > 35, then the numerator needs Int64 at least(the type we use here depends on input n)."))
    elseif typeof(n)==Int64 && n > 35
        throw(DomainError(n, "If n > 35, then the numerator needs Int128 at least (the type we use here depends on input n)."))
    elseif typeof(n)==Int128 && n > 59
        throw(DomainError(n, "If n > 59, then the numerator needs >Int128, so this code is not the code you want. Try using bernoulli(n, 0.0) to get a floating point approximation to the result."))
    elseif typeof(n)!=Int32 && typeof(n)!=Int64 && typeof(n)!=Int128 && 
        throw(DomainError(n, "n should be Int32, Int64 or Int128 at the moment"))
    end
              
    # Denominator of Bernoulli number B_n
    #   http://oeis.org/A027642
    D = [2, 6, 1, 30, 1, 42, 1, 30, 1, 66, 1, 2730, 1, 6, 1, 510, 1, 798, 1, 330, 1, 138, 1, 2730, 1, 6, 1, 870, 1, 14322, 1, 510, 1, 6, 1, 1919190, 1, 6, 1, 13530, 1, 1806, 1, 690, 1, 282, 1, 46410, 1, 66, 1, 1590, 1, 798, 1, 870, 1, 354, 1, 56786730, 1]

    # Numerator of Bernoulli number B_n (storing 62 of these because they are easy)
    #   http://oeis.org/A027641
    #    Abramowitz and Stegun, Handbook of Mathematical Function, ..., Table 23.2
    N = [-1, 1, # 2
         0, -1, # 4
         0, 1,  # 6
         0, -1, # 8
         0, 5,      # 10
         0, -691,   # 12
         0, 7,      # 14
         0, -3617,  # 16
         0, 43867,  # 18
         0, -174611,# 20
         0, 854513,        # 22    
         0, -236364091,    # 24
         0, 8553103,       # 26
         0, -23749461029,  # 28
         0, 8615841276005, # 30
         0, -7709321041217,        # 32
         0, 2577687858367,         # 34
         0, -26315271553053477373, # 36
         0, 2929993913841559,      # 38
         0, -261082718496449122051,# 40
         0, 1520097643918070802691,        # 42
         0, -27833269579301024235023,      # 44
         0, 596451111593912163277961,      # 46
         0, -5609403368997817686249127547, # 48
         0, 495057205241079648212477525,   # 50
         0, -801165718135489957347924991853,     # 52
         0, 29149963634884862421418123812691,    # 54
         0, -2479392929313226753685415739663229, # 56
         0, 84483613348880041862046775994036021, # 58
    ]
    # 0, -1215233140483755572040304994079820246041491 # 60  -- this would have to be a BigInt,
    # but most tables don't go further anyway because continuous approximations can be used
              
    # we could save some storage by removing all the 0 entries above, which aren't used, but
    # it all makes more sense to me this way 
    
    if n==0
        return Rational{typeof(n)}(1, 1)
    else
        return Rational{typeof(n)}(N[n], D[n])
    end
end


"""
    euler(n)

 Provides a lookup table for the first 60 Euler numbers ``E_n``
 e.g., see

 + Abramowitz and Stegun, Handbook of Mathematical Function, ..., Table 23.2
 + [https://en.wikipedia.org/wiki/Euler_numbers](https://en.wikipedia.org/wiki/Euler_numbers)

## Arguments
* ``n`` `::Integer`: return the nth Euler number

## Notes

The return type for this is always BigInt

## Examples
```jldoctest; setup = :(using Polylogarithms)
julia> bernoulli(6)
1//42
```

"""
function euler(n::Integer)

    if n<0
        throw(DomainError(n))
    elseif n > 1 && isodd(n)
        return BigInt(0)
    elseif n > 60
        throw(DomainError(n, "Largest output value here is n=60."))
    end

    E = [1, # 0
         0, -1, 
         0, 5,
         0, -61, 
         0, 1385, 
         0, -50521, # 10
         0, 2702765,
         0, -199360981,
         0, 19391512145,
         0, -2404879675441,
         0, 370371188237525, # 20
         0, -69348874393137901,
         0, 15514534163557086905,
         0, -4087072509293123892361,
         0, 1252259641403629865468285,
         0, -441543893249023104553682821, # 30
         0, 177519391579539289436664789665,
         0, -80723299235887898062168247453281,
         0, 41222060339517702122347079671259045,
         0, -23489580527043108252017828576198947741,
         0, 14851150718114980017877156781405826684425, # 40
         0, -10364622733519612119397957304745185976310201,
         0, 7947579422597592703608040510088070619519273805,
         0, -6667537516685544977435028474773748197524107684661,
         0, 6096278645568542158691685742876843153976539044435185,
         0, -6053285248188621896314383785111649088103498225146815121, # 50
         0, 6506162486684608847715870634080822983483644236765385576565,
         0, -7546659939008739098061432565889736744212240024711699858645581,
         0, 9420321896420241204202286237690583227209388852599646009394905945,
         0, -12622019251806218719903409237287489255482341061191825594069964920041,  
         0, 18108911496579230496545807741652158688733487349236314106008095454231325 # 60
    ]
    # we could save some storage by removing all the 0 entries above, which aren't used, but
    # it all makes more sense to me this way 

    return E[n+1]
end
