# test using data
include("test_defs.jl")

L = Symbol("Li_s(z)")
@testset "Polylogarithm on random test data (validated using Mathematica)" begin
    for C=1:4
        filename = @sprintf("../data/polylog_test_data_rand_%d.csv", C)
        data1 = CSV.read(filename, DataFrame; delim=",", types=String)

        # has trouble reading in numbers like "2." so read all into strings, and parse
        data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
        data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
        data1[!, L] = parse.(Complex{Float64}, data1[!, L] )
        
        m = size(data1,1)
        Li = data1[!, L]
        s  = data1[!,:s]
        z  = data1[!,:z]
        no_failed = 0
        
        for i=1:m
            rel_error =  ( polylog(s[i],z[i]) - Li[i] )./ Li[i]
            if abs(rel_error) > accuracy_goal1
                no_failed += 1
                print("   accuracy warning: C=$C, s=$(s[i]), z=$(z[i]), relative error = ")
                println( abs(rel_error) )
            end
        end
        @test no_failed < 10 # NB don't expect all these to pass, just most
    end
end

@testset "Dilogarithm on test data from Morris" begin
    # note that Morris computes only the real values

    # 100 tests

    filename = "../data/morris_dilog_appendixB.csv"
    data2 = CSV.read(filename, DataFrame; delim=",", comment="#", types=Float64)
    
    m = size(data2,1) 
    Lx = Symbol("dilog(x)")
    Lx_m = Symbol("dilog(-x)")

    x = data2[:, :x]
    Li = polylog.(2, x)
    Li_m = polylog.(2, .- x)

    branch = zeros(m)
    k = findall(x .>= 1.0  )
    branch[k] = -pi .* log.(x[k]) .^ 1  ./ gamma(2)

    error1 = Li .-  ( data2[:, Lx] .+ branch * im )
    error2 = Li_m .-  data2[:, Lx_m]
    rel_error1 =  abs.( error1 ./ data2[:, Lx] )
    rel_error2 =  abs.( error2 ./ data2[:, Lx_m] )

    max_rel_error1 = maximum( rel_error1 )
    max_rel_error2 = maximum( rel_error2 )

    no_failed1 = sum( rel_error1 .> accuracy_goal1 )
    no_failed2 = sum( rel_error2 .> accuracy_goal1 )

    # k1 = findall( rel_error1 .> accuracy_goal1 )
    # k2 = findall( rel_error2 .> accuracy_goal1 )
    # for i=k1
    #     print("   error warning: s=2, z=$(x[i]), relative error = ")
    #     println( abs(rel_error1[i]) )
    # end
    # for i=k2
    #     print("   error warning: s=2, z=$(-x[i]), relative error = ")
    #     println( abs(rel_error2[i]) )
    # end
    # @test no_failed1 <= 0
    # @test no_failed2 <= 0

    for i=1:m
        @test rel_error1[i] <= accuracy_goal1
        @test rel_error2[i] <= accuracy_goal1
    end
end

