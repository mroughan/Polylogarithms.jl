"""
    rogers(z)

 Calculates Rogers L-function
     https://mathworld.wolfram.com/RogersL-Function.html
    Directly related to the dilogartihm.

 Note there are two possible versions, we use the one Bytsko, 1999, https://arxiv.org/abs/math-ph/9911012
 but note that this is a normalised form (the extra (6/π^2) such that rogers(1)=1) as compared to Rogers (1907)
"""
function rogers(z::Number;
                 level=1, # keep track of recursion
                 accuracy::Float64=default_accuracy,
                 min_iterations::Integer=0,
                max_iterations::Integer=default_max_iterations)
    Li = polylog(2.0, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
    small_number = 1.0e-12
    if abs(z) < small_number
        return 0.0 + 0.0im
    elseif abs(1-z) < small_number
        return 1.0 + 0.0im
    else
        z = convert(Complex{Float64}, z)
        return (6/π^2) * ( Li + log(z)*log(1-z)/ 2.0 )
    end
end

"""
    spence(z)

    An alias for the dilogarith (polylog with s=2)
"""
function spence(z::Number;
                level=1, # keep track of recursion
                accuracy::Float64=default_accuracy,
                min_iterations::Integer=0,
                max_iterations::Integer=default_max_iterations)
    return polylog(2.0, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
end


"""
    dilog(z)

    Just alias for the polylogarith with s=2
"""
function dilog(z::Number;
                level=1, # keep track of recursion
                accuracy::Float64=default_accuracy,
                min_iterations::Integer=0,
                max_iterations::Integer=default_max_iterations)
    return polylog(2.0, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
end

"""
    trilog(z)

    Just alias for the polylogarith with s=3
"""
function trilog(z::Number;
                level=1, # keep track of recursion
                accuracy::Float64=default_accuracy,
                min_iterations::Integer=0,
                max_iterations::Integer=default_max_iterations)
    return polylog(3.0, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
end

"""
    tetralog(z)

    Just alias for the polylogarith with s=4
"""
function tetralog(z::Number;
                level=1, # keep track of recursion
                accuracy::Float64=default_accuracy,
                min_iterations::Integer=0,
                max_iterations::Integer=default_max_iterations)
    return polylog(4.0, z; level=level, accuracy=accuracy, min_iterations=min_iterations, max_iterations=max_iterations)
end
