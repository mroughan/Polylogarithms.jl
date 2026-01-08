# test using data
include("test_defs.jl")

Z = Symbol("zeta'(s)")
@testset "Zeta derivatives on random test data (validated using Mathematica)" begin
    filename = @sprintf("../data/zeta_derivative_data_rand.csv")
    data1 = CSV.read(filename, DataFrame; delim=",", types=String)
    
    # has trouble reading in numbers like "2." so read all into strings, and parse
    data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
    data1[!, Z] = parse.(Complex{Float64}, data1[!, Z] )
    
    m = size(data1,1)
    ζ_d = data1[!, Z]
    s  = data1[!,:s]
    no_failed = 0
    rel_error = zeros(Complex{Float64}, m)

    for i=1:m
        rel_error[i] =  ( zeta_derivative(s[i])[1] - ζ_d[i] )./ ζ_d[i]
        if abs(rel_error[i]) > accuracy_goal1
            no_failed += 1
            print("   accuracy warning: s=$(s[i]), relative error = ")
            println( abs(rel_error[i]) )
        end
    end
    @test no_failed <= 0 # NB don't expect all these to pass, just most
    max_error = maximum(abs.(rel_error))
    @test max_error < 1.0e-11
    println("   maximum error = $max_error  and no_failed = $no_failed ")
end

