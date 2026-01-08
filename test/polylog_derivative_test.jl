# test using data
include("test_defs.jl")

L = Symbol("dLi_s/dz")
@testset "Polylogarithm derivative with respect to z on random test data (validated using Mathematica)" begin
    filename = @sprintf("../data/polylog_test_data_rand_dz.csv")
    data1 = CSV.read(filename, DataFrame; delim=",", types=String)

    # has trouble reading in numbers like "2." so read all into strings, and parse
    data1[!,:s] = parse.(Complex{Float64}, data1[!,:s] )
    data1[!,:z] = parse.(Complex{Float64}, data1[!,:z] )
    data1[!, L] = parse.(Complex{Float64}, data1[!, L] )
    
    m = size(data1,1)
    dLidz = data1[!, L]
    s  = data1[!,:s]
    z  = data1[!,:z]
    no_failed1 = 0
    rel_error1 = zeros(Complex{Float64}, m)
    
    for i=1:m
        rel_error1[i] =  ( polylog_dz(s[i],z[i]) - dLidz[i] )./ dLidz[i]
        if abs(rel_error1[i]) > accuracy_goal1
            no_failed1 += 1
            print("   accuracy warning: s=$(s[i]), z=$(z[i]), relative error = ")
            println( abs(rel_error1[i]) )
        end
    end
    @test no_failed1 < 10 # NB don't expect all these to pass, just most
end 

L = Symbol("dLi_s/ds")
@testset "Polylogarithm derivative with respect to s on random test data (validated using Mathematica)" begin
    filename = @sprintf("../data/polylog_test_data_rand_ds.csv")
    data2 = CSV.read(filename, DataFrame; delim=",", types=String)

    # has trouble reading in numbers like "2." so read all into strings, and parse
    data2[!,:s] = parse.(Complex{Float64}, data2[!,:s] )
    data2[!,:z] = parse.(Complex{Float64}, data2[!,:z] )
    data2[!, L] = parse.(Complex{Float64}, data2[!, L] )
    
    m = size(data2,1)
    dLids = data2[!, L]
    s  = data2[!,:s] 
    z  = data2[!,:z]
    no_failed2 = 0
    rel_error2 = zeros(Complex{Float64}, m)
    
    for i=1:m
        output = polylog_ds(s[i],z[i])
        rel_error2[i] =  (output[1]  - dLids[i] )./ dLids[i]
        if abs(rel_error2[i]) > accuracy_goal1
            no_failed2 += 1
            print("   accuracy warning: s=$(s[i]), z=$(z[i]), k=$(output[2]), series=$(output[3]), relative error = ")
            println( abs(rel_error2[i]) )
        end 
    end
    @test no_failed2 < 180 # NB don't expect all these to pass, just most, but derivatives are harder
    max_error2 = maximum(abs.(rel_error2))
    @test max_error2 < 1.0e-9
    println("   maximum error = $max_error2  and no_failed = $no_failed2 ")

# scatter( -real.(s), log10.( abs.(rel_error2) ) ; label="", xlabel="-real(s)", ylabel="log10(abs(rel error))", ylim = (-13, -9), markersize =4, markeralpha=0.3)

# scatter( abs.(z), log10.( abs.(rel_error2) ) ; label="", xlabel="-real(s)", ylabel="log10(abs(rel error))", ylim = (-13, -9), markersize =4, markeralpha=0.3)

end

0
