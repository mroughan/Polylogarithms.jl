include("test_defs.jl")

@testset "Harmonic numbers" begin

    @testset "    throws errors" begin
        @test_throws DomainError harmonic(-1)
        @test_throws MethodError harmonic(Float32(1.0))
    end
    
    @testset "    types" begin
        @test typeof(harmonic(1)) == Float64
        @test typeof(harmonic(1.0)) == Float64
        @test typeof(harmonic(complex(1.0))) == Complex{Float64}
    end

    @testset "    basics" begin
        @test harmonic(1) ≅ 1.0
        @test harmonic(2) ≅ 1.5
        @test harmonic(3) ≅ 11.0/6.0
        @test harmonic(4) ≅ 25.0/12.0
        @test harmonic(5) ≅ 137.0/60.0
    end

    @testset "    identities" begin
        for n=4:10
            @test harmonic(n) ≅ 1.0/n + harmonic(n-1)
        end
        
        for n=20:30
            @test harmonic(n) ≅ sum( 1.0./(1:n) )
        end
    end

    @testset "    dataset 1 (integer arguments)" begin
        data1 = CSV.read("../data/harmonic_test_data_1.csv", DataFrame)
        m = size(data1,1)
        H = Symbol("real(H_n)")
        for i=1:m
            @test harmonic( data1[i,:n] ) ≅  data1[i,H]
        end
    end
 
    @testset "    dataset 2 (real arguments)" begin
        data2 = CSV.read("../data/harmonic_test_data_2.csv", DataFrame)
        m = size(data2,1)
        H = Symbol("real(H_n)")
        for i=1:m
            @test harmonic( data2[i,:x] ) ≈  data2[i,H] # note the lower accuracy requirement here
        end
    end

    @testset "    dataset 3 (complex arguments)" begin
        data3 = CSV.read("../data/harmonic_test_data_3.csv", DataFrame)
        m = size(data3,1)
        zr = Symbol("Re(z)")
        zi = Symbol("Im(z)")
        Hr = Symbol("Re(H(z))")
        Hi = Symbol("Im(H(z))")
        for i=1:m
            z = data3[i,zr] + im*data3[i,zi]
            Hz = data3[i,Hr] + im*data3[i,Hi]
            @test harmonic( z ) ≈  Hz # note the lower accuracy requirement here
        end
    end
end

@testset "Generalised Harmonic numbers" begin
      
    @testset "    throws errors" begin
        @test_throws DomainError harmonic(-1, 1.2)
        @test_throws MethodError harmonic(1.0, 1)
        @test_throws MethodError harmonic(1.0, 1.0)
        @test_throws MethodError harmonic(complex(1.0), 1)
    end

    @testset "    output types" begin
        @test typeof(harmonic(1, 1)) == Float64
        @test typeof(harmonic(1, 1.0)) == Float64
    end

    @testset "     identities" begin
        for n=1:10
            @test harmonic(n,1.0) ≅ harmonic(n)
            @test harmonic(n,0.0) ≅ n
            
            for r = 3:2:9
                @test harmonic(n,r) ≅ Float64(n)^(-r) + polygamma(r-1,n)/gamma(r) + zeta(r)
            end
        end
    end

    @testset "    dataset 4 (integer arguments)" begin
        H = Symbol("real(H_{n,r})")
        data4 = CSV.read("../data/harmonic_test_data_4.csv", DataFrame)
        m = size(data4,1)
        for i=1:m
            @test harmonic( data4[i,:n], data4[i,:r] ) ≅ data4[i,H]
        end
    end

    @testset "    dataset 5 (integer n, real r)" begin
        H = Symbol("real(H_{n,r})")
        data5 = CSV.read("../data/harmonic_test_data_5.csv", DataFrame)
        m = size(data5,1)
        for i=1:m
            @test harmonic( data5[i,:n], data5[i,:r] )  ≅ data5[i,H]
        end
    end
end
