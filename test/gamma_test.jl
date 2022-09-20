include("test_defs.jl")
include("../src/gamma_derivatives.jl")
# plan to do more here

@testset "Derivatives of the gamma function at 1.0" begin

    @testset "    throws errors" begin
        @test_throws DomainError D(2, -1, 1.0)
        @test_throws DomainError D(-1, 2, 1.0)
        # @test_throws MethodError stieltjes(1.5)
    end
    
    @testset "    types" begin
        # @test typeof(stieltjes(1)) == Float64
    end

    @testset "    values" begin
        @test g1(0) ≈ 1.0
        @test g2(0) ≈ 1.0
        @test g1(1) ≈ -γ
        @test g2(1) ≈ -γ
        for n=2:5
            @test g1(n) ≈ g2(n)
        end
    end
end
