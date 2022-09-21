include("test_defs.jl")

@testset "Stieltjes numbers" begin

    @testset "    throws errors" begin
        @test_throws DomainError stieltjes(-1)
        @test_throws DomainError stieltjes(79)
        @test_throws MethodError stieltjes(1.5)
    end
    
    @testset "    types" begin
        @test typeof(stieltjes(1)) == Float64
    end

    @testset "    values" begin
        @test stieltjes(0) ≅ γ
        function gen_euler_calc(n::Integer)
            m = 10000000 # just choose a big value for the moment
            total = -log(m)^(n+1) / (n+1)
            for k=1:m
                total += log(k)^n / k
            end
            return total
        end
        @test abs(stieltjes(1) - gen_euler_calc(1)) < 1.0e-5 # this series is a little crude as a test
        @test abs(stieltjes(2) - gen_euler_calc(2)) < 1.0e-4 # this series is a little crude as a test
        @test abs(stieltjes(3) - gen_euler_calc(3)) < 1.0e-3 # this series is a little crude as a test
    end

   
    @testset "    values (2)" begin
        # these values are taken from alt sources to provide some tests
        
        #   https://beta.lmfdb.org/riemann/stieltjes/
        values1 = CSV.read("../data/stieltjes_lmfdb.csv", DataFrame; comment="#")

        #   https://fredrikj.net/math/stieltjes100k20d.txt
        values2 = CSV.read("../data/stieltjes100k20d.txt", DataFrame; comment="#")
                  
        for i=0:78
            s = stieltjes(i)
            if i != 55
                # 55th element of the LMFDB data lacks a digit of precision
                @test s == values1[i+1,:gamma_k]
            end
            @test s == values2[i+1,:gamma_k]
        end
    end
end

