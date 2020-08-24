# test functionality and identities
using Polylogarithms
using SpecialFunctions
using Test
using DataFrames, CSV
import Base.MathConstants: π, pi, ℯ, e, γ, eulergamma, catalan, φ, golden
include("test_defs.jl")
include("../bench/utilities.jl")
Q = Polylogarithms.Q

@testset "Polylogarithm polylog" begin
    
    # check throws errors
    @testset "    throws errors" begin
        # @test_throws DomainError polylog(1, 0)              # should work for all inputs
        # @test_throws MethodError polylog( 1, Float32(0.0))  # should work for all numbers
    end 

    # check output types
    @testset "    output types" begin
        @test typeof( polylog(1, 0) ) == Float64
        @test typeof( polylog(complex(1), 0) ) == Complex{Float64}
        # we need to do some more testing here and consider carefully how consistent this should be
    end
 
    # check it can handle all types of inputs
    @testset "    input formats" begin
        @test polylog(-1.0, 0.0) ≈ 0.0
        @test polylog(-1,   0.0) ≈ 0.0
        @test polylog(-1.0, Complex(0.0)) ≈ 0.0
        @test polylog(Complex(-1.0), 0.0) ≈ 0.0
        @test polylog(Complex(-1.0), 0.0) ≈ 0.0
        @test polylog(1, 0) ≈ 0.0
        @test polylog(1.0, 0) ≈ 0.0
        @test polylog(1, 0.0) ≈ 0.0
        @test polylog(1.0, 0.0) ≈ 0.0
        @test polylog(1 // 2, 0) ≈ 0.0
        @test polylog(1 // 2, 1 // 2) ≈ 0.8061267230428526
        @test polylog( π, π ) ≈ 3.920844767506471 - 1.8340224080910017im
        x = collect(0.0:0.1:0.9)
        @test all([polylog.(1, x)[i] ≈ polylog(1, x[i]) for i=1:length(x)])
        @test polylog(Complex(-1.0), Complex(0.3)) ≈ polylog(-1.0, 0.3)
        @test polylog(Complex(-1.0), Complex(0.3)) ≈ polylog(-1.0, 0.3)
    end

    # check specific known values
    @testset "     s = n (a real integer)" begin
        # simple cases
        @test polylog(1, 0.5) ≈ log(2)
        @test !isfinite( polylog(0.99999, 1.0) ) # the singularity
        @test !isfinite( polylog(-1, 1.0) )      # the singularity
        @test !isfinite( polylog(0.5 + 0.5im, 1.0) )     # the singularity
        @test polylog(1, 2) ≈ -pi*im # comes from -log(complex(1-2))
    
        @testset "     dilogarithm for real z" begin
            @test polylog(2,-1.0)    ≈ -pi^2/12.0
            @test polylog(2, 0.0)    ≈ 0.0
            @test polylog(2, 0.5)    ≈ pi^2/12 - 0.5*log(2)^2
            @test polylog(2, 1.0)    ≈ pi^2/6.0
            @test polylog(2, 2.0)    ≈ pi^2/4 - im*pi*log(2)
            @test polylog(2, -φ)     ≈ -pi^2/10 - log(φ)^2
            @test polylog(2, -1/φ)   ≈ -pi^2/15 + log(φ)^2/2
            @test polylog(2, 1/φ^2)  ≈  pi^2/15 - log(φ)^2
            @test polylog(2, 1/φ)    ≈  pi^2/10 - log(φ)^2
            @test polylog(2, φ)      ≈  11*pi^2/15 + log(Complex(-1/φ))^2/2 # wiki has this one, but no ref
            @test polylog(2, φ^2)    ≈ -11*pi^2/15 - log(Complex(-φ))^2

            # identities for dilogarithm
            Z = [3.0 + 0.4im, -3.0 + 0.4im, 3.0 - 0.4im, -3.0 + -0.4im]
            for i=1:length(Z)
                z = Z[i]
                @test polylog(2, z) + polylog(2, 1/z) ≈ -pi^2/6.0 - log(Complex(-z))^2/2.0
            end
        end

        @testset "     trilogarithm for real z" begin
            # https://mathworld.wolfram.com/Trilogarithm.html
            @test polylog(3,-1.0)             ≈ -3*zeta(3)/4
            @test polylog(3, 0.0)             ≈ 0.0
            @test polylog(3, 0.5)             ≈ log(2)^3/6.0 - pi^2*log(2)/12.0 + (7.0/8.0)*zeta(3)
            @test polylog(3, 1.0)             ≈ zeta(3)
            @test polylog(3, Float64(φ)^(-2)) ≈ 4*zeta(3)/5 + 2*log(φ)^3/3 - 2*pi^2*log(φ)/15
        end
        
        @testset "     general case for real z" begin
            X = collect(-3.0:0.1:3.0)
            for i=1:length(X)
                x = X[i]
                # println("x = $x")
                @test polylog(1, x) ≈ -log(Complex(1-x))
                @test polylog(0, x) ≈ x ./ (1-x)
                @test polylog(-1, x) ≈ x ./ (1-x).^2
                @test polylog(-2, x) ≈ x .* (1+x) ./ (1-x).^3
                @test polylog(-3, x) ≈ x .* (1+4*x+x.^2) ./ (1-x).^4
                @test polylog(-4, x) ≈ x .* (1+x) .* (1+10*x+x.^2) ./ (1-x).^5
            end
        end
        
        @testset "     general case for complex z" begin
            X = collect(-3.0:0.5:3.0)
            Y = [-1.3, -0.4, 0.4, 1.5]
            for i=1:length(X)
                for j=1:length(Y)
                    z = Complex(X[i], Y[j])
                    # println("z = $z")
                    @test polylog(1, z) ≈ -log(Complex(1-z))
                    @test polylog(0, z) ≈ z ./ (1-z)
                    @test polylog(-1, z) ≈ z ./ (1-z).^2
                    @test polylog(-2, z) ≈ z .* (1+z) ./ (1-z).^3
                    @test polylog(-3, z) ≈ z .* (1+4*z+z.^2) ./ (1-z).^4
                    @test polylog(-4, z) ≈ z .* (1+z) .* (1+10*z+z.^2) ./ (1-z).^5
                end
            end
        end    
    end

    @testset "    particular values |z| == 1" begin
        S_r = [2.1 2.5 3.0]
        S_i = [-1.3, -1.0, -0.5, 0.0, 0.5, 1.0, 1.3]
        for i=1:length(S_r)
            for j=1:length(S_i)
                s = Complex(S_r[i], S_i[j])
                @test polylog(s, -1.0) ≈ -eta(s)
                @test polylog(s,  1.0) ≈  zeta(s)
                @test polylog(s,  im) ≈  - 2.0^(-s)*eta(s) + im*dirichlet_beta(s) 
            end
        end
    end

    @testset "    additional Identities" begin
        z = 0.5
        for n=1:5
            @test polylog(-n,z) + (-1)^n * polylog(-n, 1/z) ≈ 0.0 
        end

        # for real s, and real z<1, polylog should be real
        S = [-1, 0.1, 2]
        Z = [-2, -1.0, 0.1, 0.95]
        for i=1:length(S)
            for j=1:length(Z)
                s = S[i]
                z = Z[j]
                # println("s = $s; z = $z")
                @test abs( imag( polylog(s,  z) ) ) < 1.0e-12 # actually they usually do better than this
            end
        end

        # for real s, and real z>=1, the imaginary part is given
        S = [-1.5, 0.1, 2]
        Z = [1.05, 3.0]
        for i=1:length(S)
            for j=1:length(Z)
                s = S[i]
                z = Z[j]
                # println("s = $s; z = $z")
                μ = log(z)
                @test imag( polylog(s,  z) ) ≈ -pi*μ^(s-1)/gamma(s)
            end
        end

        # duplication formula (or square formula)
        S_r = [2.1 2.5 3.0]
        S_i = [-1.3, -1.0, -0.5, 0.0, 0.5, 1.0, 1.3]
        z = 3.0 - 0.1im
        for i=1:length(S_r)
            for j=1:length(S_i)
                s = Complex(S_r[i], S_i[j])
                # println("s = $s; z = $z")
                @test polylog(s,z) + polylog(s,-z) ≈ complex(2)^(1-s) * polylog(s, z^2)
            end
        end
    end
end
