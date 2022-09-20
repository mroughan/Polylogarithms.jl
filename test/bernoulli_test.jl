include("test_defs.jl")

@testset "Bernoulli numbers" begin
    @testset "    throws errors" begin
        @test_throws DomainError bernoulli(-1)
        @test_throws DomainError bernoulli( Int128(60) )
        @test_throws DomainError bernoulli( Int64(36)  )
        @test_throws DomainError bernoulli( Int32(24)  )
        @test_throws DomainError bernoulli( BigInt(1)  )
        @test_throws DomainError bernoulli( Int8(1)  )
        @test_throws MethodError bernoulli(1.0)
    end

    @testset "    types" begin
        @test typeof( bernoulli( Int32(2)  )) == Rational{Int32}
        @test typeof( bernoulli( Int64(2)  )) == Rational{Int64}
        @test typeof( bernoulli( Int128(2)  )) == Rational{Int128}
    end
    
    @testset "    approximation for explicit cases" begin
        for i=2:34
            @test Float64( bernoulli(i) ) ≈ -i*zeta( 1 - i )
        end
        for i=35:59
            @test Float64( bernoulli( Int128(i) ) ) ≈ -i*zeta( 1 - i )
        end
    end

    @testset "    dataset 1" begin
        B = Symbol("B_n")
        data1 = CSV.read("../data/bernoulli_test_data_1.csv", DataFrame)
        m = size(data1,1)
        for i=1:m
            @test bernoulli( data1[i,:n] ) ≅ parse(Rational{Int64},data1[i,B])
        end
    end

    @testset "    identities" begin
        # sum formula
        for n=2:35
            k = collect( 0:Int64(n) )
            @test sum( binomial.(n+1, k) .*  bernoulli.(k) ) == 0
        end
    end
end

@testset "Euler numbers" begin
    @testset "    throws errors" begin
        @test_throws DomainError euler(-1)
        @test_throws DomainError euler( 62 )
    end

    @testset "    types" begin
        @test typeof( euler(2) ) == BigInt
    end
    
    @testset "    congruences" begin
        for p in [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]
            if mod(p,4) == 1
                @test mod( (-1)^( (p-1)/2 ) * euler(p-1), p) == 0
            elseif mod(p,3) == 3
                @test mod( (-1)^( (p-1)/2 ) * euler(p-1), p) == p-2
            else
                # ignore these cases
            end
        end
    end
    
end

@testset "Bernoulli polynomials" begin
    @testset "    throws errors" begin
        @test_throws DomainError bernoulli(-1, 0.0)
        @test_throws MethodError bernoulli(1.0, 0.0)
        @test_throws MethodError bernoulli(1.0, 1)
    end

    @testset "    types" begin
        @test typeof(bernoulli(2, 1.0)) == Float64
    end
 
    @testset "    consistency with Bernoulli numbers" begin
        for i=0:14
            @test Float64(bernoulli(i)) ≅ bernoulli(i, 0.0)
        end
        for i=15:34
            # somewhat reliant here on accuracy of SpecialFunctions, so can't make these better without rewriting that
            @test Float64(bernoulli(i)) ≈ bernoulli(i, 0.0)
        end
     end
    
    @testset "    explicit cases" begin
        X = collect(-2.0 : 0.1 : 3.0)
        for i=1:length(X)
            x = X[i]
            # println("   x = $x")
            @test bernoulli(0, x) ≅ 1.0
            @test bernoulli(1, x) ≅ x   - 0.5
            @test bernoulli(2, x) ≅ x^2 -     x   + 1.0/6.0
            @test bernoulli(3, x) ≅ x^3 - 1.5*x^2 + 0.5*x
            @test bernoulli(4, x) ≅ x^4 - 2.0*x^3 +     x^2 - 1/30.0
            @test bernoulli(5, x) ≅ x^5 - 2.5*x^4 + (5.0/3.0)*x^3 - x/6.0
            # @test bernoulli(6, x) ≈ x^6 - 3.0*x^5 + (5.0/2.0)*x^4 - 0.5*x^2 + 1/42.0
            # Hurwitz-zeta doesn't work well enough here
        end
    end
    
    @testset "    symmetry identities" begin
        X = collect(0.0: 0.1 : 1.0)
        for i=1:length(X)
            x = X[i]
            # for n=1:9
            for n=1:6
                # println("   n = $n, x = $x")
                # start to get errors of order 1.0e-14 for n=5, 1.0e-13 around n=10, ...
                @test bernoulli(n, 1-x)  ≈  (-1.0)^n * bernoulli(n, x)
                @test bernoulli(n, x+1)  ≈  bernoulli(n, x) + n*x^(n-1)
                # SpecialFunctions.zeta(-6, -1.0) = NaN (probably should be 2)
                @test (-1)^n * bernoulli(n, -x)  ≈   bernoulli(n, x) + n*x^(n-1) 
                
                # Raabe (1851)
                m = 6.0
                k = 0:m-1
                @test sum( bernoulli.(n, x .+ k./m ) )/m  ≈  m^(-n)*bernoulli(n, m*x )
            end
        end
    end
 
    @testset " dataset 2" begin 
        B = Symbol("B_n(x)")
        data2 = CSV.read("../data/bernoulli_test_data_2.csv", DataFrame; delim=",", types=String)
        data2[!,:n] = parse.(Float64, data2[!,:n] )
        data2[!,:x] = parse.(Float64, data2[!,:x] )
        data2[!,B] = parse.(Float64, data2[!,B] )
        m = size(data2,1)
        for i=1:m
            @test bernoulli( Int(data2[i,:n]),  data2[i,:x]) ≈ data2[i,B]
        end
    end
    
    @testset " dataset 3" begin
        B = Symbol("B_n(x)")
        data3 = CSV.read("../data/bernoulli_test_data_3.csv", DataFrame; delim=",")
        m = size(data3,1)
        for i=1:m
            @test bernoulli( Int(data3[i,:n]),  data3[i,:x]) ≈ data3[i,B]
        end
   end    
end

