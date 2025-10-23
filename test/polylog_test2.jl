# test using data
include("test_defs.jl")

L = Symbol("Li_s(z)")
@testset "Polylogarithm on test data sets -- NB don't expect all these to pass, just most" begin
    for C=1:3
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

        for i=1:m
            rel_error =  ( polylog(s[i],z[i]) - Li[i] )./ Li[i]
            if abs(rel_error) > accuracy_goal1
                print("   error warning: C=$C, s=$(s[i]), z=$(z[i]), relative error = ")
                println( abs(rel_error) )
            end
        end

    end
end
