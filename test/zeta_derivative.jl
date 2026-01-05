include("test_defs.jl")
include("../src/zeta_derivative.jl")

@testset "Riemann zeta derivatives" begin

    @testset "    throws errors" begin
        @test_throws DomainError stieltjes(-1)
        @test_throws DomainError stieltjes(79)
        @test_throws MethodError stieltjes(1.5)
    end
    
    @testset "    types" begin
        @test typeof(stieltjes(1)) == Float64
    end

    @testset "    values" begin
        # https://mathworld.wolfram.com/RiemannZetaFunction.html
        @test zeta_derivative(2) ≅ (1/6) * π^2 * ( γ + log(2*π) - 12*log(A) )
                                  # = -0.93754825431...
        
        @test zeta_derivative(0.5) ≅ (1/4) * ( π + 2*γ + 6*log(2) + 2*log(π) ) * zeta(0.5)
                                  # = -3.92264613...

        @test zeta_derivative(-2.0) ≅ -zeta(3) / (4*π^2)
        @test zeta_derivative(-4.0) ≅ 3 * zeta(5) / (4*π^4)
        @test zeta_derivative(-6.0) ≅ -45 * zeta(7) / (8*π^6)
        @test zeta_derivative(-8.0) ≅ 315 * zeta(9) / (4*π^8)
        # could use general form for ζ(-2n), see https://mathworld.wolfram.com/RiemannZetaFunction.html
        
        @test zeta_derivative(-1) ≅ (1/12) - log(A)

        @test zeta_derivative(0.0) ≅ -(1/2) * log(2 π)
 
        
        
    end

end
 
