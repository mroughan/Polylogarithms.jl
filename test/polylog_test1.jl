# test functionality and identities
include("test_defs.jl")
Q = Polylogarithms.Q

@testset "Polylogarithm polylog identities and special values" begin
    
    # check throws errors
    @testset "    throws errors" begin
        # @test_throws DomainError polylog(1, 0)              # should work for all inputs
        # @test_throws MethodError polylog( 1, Float32(0.0))  # should work for all numbers
        @test_throws MethodError polylog( 1, "1.0")  # should work for all numbers
        @test_throws MethodError polylog( "1.0", 1)  # should work for all numbers
        @test_throws MethodError polylog( 1 )
        # @test_throws ArgumentError polylog( 1 )
    end 

    # check output types
    @testset "    output types" begin
        @test typeof( polylog(1, 0) ) == Float64
        @test typeof( polylog(complex(1), 0) ) == Complex{Float64}
        # we need to do some more testing here and consider carefully how consistent this should be
    end
 
    # check keywords work
    @testset "    output types" begin
        @test polylog(-1, 0; level=1, accuracy=1.0e-5,  min_iterations=1, max_iterations=100) ≈ 0.0
    end
 
    # check subsidary routine errors, just to make coverage cleaner
    @testset "    unexported function errors" begin
        @test_throws DomainError Polylogarithms.polylog_series_1(1.0, 2.0)
        @test_throws DomainError Polylogarithms.polylog_series_1(1.0, 0.75)
        @test_throws DomainError Polylogarithms.polylog_series_2(1.0, 0.0)
        @test_throws DomainError Polylogarithms.polylog_series_3(1.0, 0.0)
        @test_throws DomainError Polylogarithms.polylog_series_3(-1.0, 0.6)
        @test_throws DomainError Polylogarithms.g_crandall(-1)
        @test_throws DomainError Polylogarithms.g_crandall(15)
        τ = 0.00001
        @test_throws DomainError Polylogarithms.Q_closed(0, τ, 0; n_terms=0)
        @test_throws DomainError Polylogarithms.Q_closed(0, τ, 0; n_terms=4)
        @test_throws DomainError Polylogarithms.Q(0, τ, 0; n_terms=0)
        @test_throws DomainError Polylogarithms.Q(0, τ, 0; n_terms=8)
     end

    @testset "    unexported function values" begin
        ℒ = 0.1
        d2 = SpecialFunctions.digamma(1) - ℒ
        @test Polylogarithms.c_closed(0, 0, ℒ) ≈ harmonic(0) - ℒ
        @test Polylogarithms.c_closed(0, 1, ℒ) ≈ -stieltjes(1)   - d2^2/2      - π^2/6 + SpecialFunctions.polygamma(1,1)/2
        @test Polylogarithms.c_closed(0, 2, ℒ) ≈  stieltjes(2)/2 + d2^3/6 + d2*( π^2/6 -
                                                       SpecialFunctions.polygamma(1,1)/2 ) +
                                                       SpecialFunctions.polygamma(2,1)/6
        for n=0:3
            @test Polylogarithms.c_closed(n, 0, ℒ) ≈ Polylogarithms.c_crandall(n, 0, ℒ)
            @test Polylogarithms.c_closed(n, 1, ℒ) ≈ Polylogarithms.c_crandall(n, 1, ℒ)
            @test Polylogarithms.c_closed(n, 2, ℒ) ≈ Polylogarithms.c_crandall(n, 2, ℒ)
        end
        for k = 1:5
            @test Polylogarithms.f_crandall(k,1) == -harmonic(k)
        end
        @test Polylogarithms.g_crandall(1) ≈ -γ

        τ = 0.00001
        n = 1
        @test Polylogarithms.Q_closed(n, τ, ℒ; n_terms=3) == Polylogarithms.Q(n, τ, ℒ; n_terms=3)
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

    # cases that caused errors in the past
    @testset "     cases to check old errors don't recur" begin
        @test polylog( 0.000964487315968654, 0.25 ) ≈  0.3332677049928309 # small, positive s shouldn't use series 3
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
    
    # diagnostics
    @testset "    diagnostics" begin
        S = [ -1.0 ,  0.5+0.5im,  3.0,   -1.5,         -1.5,        -1.5,      1   ]
        Z = [ -0.25, 1.0 + im,   -1.0,  -50.0, -400 + 250im, 600 + 600im,  -500.0  ]

        # series 1 cases, polylog(-1, z) ≈ z ./ (1-z).^2
        i = 1
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( z ./ (1-z).^2, 24, 1, 0) )

        # series 2 cases
        i = 2
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 16, 2, 0) )
        
        # series 3 cases,    polylog(3,-1.0)             ≈ -3*zeta(3)/4
        i = 3
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( -3*zeta(3)/4, 36, 3, 0) )
        
        # m-th root multiplication formula
        i = 4
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 36, 22, 0) )
      
        i = 5
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 58, 32, 0) )

        i = 6
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z),  52, 42, 0) )

        # duplication recursion
        i = 7
        s = S[i]
        z = Z[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 120, ((3, 3), (3, 3)), 1) )
    end
                   
end
