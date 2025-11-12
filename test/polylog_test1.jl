# test functionality and identities
include("test_defs.jl")
Q = Polylogarithms.Q
Z = [3.0 + 0.4im, -3.0 + 0.4im, 3.0 - 0.4im, -3.0 + -0.4im, 2.0 + 0.1im, -5.0 + 0.1im, 0.0 + 6im]

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

    # check specific known values, and consistency of identities
    @testset "     s = n (a real integer)" begin
        # simple cases
        @test polylog(1, 0.5) ≈ log(2)
        @test !isfinite( polylog(0.99999, 1.0) ) # the singularity
        @test !isfinite( polylog(-1, 1.0) )      # the singularity
        @test !isfinite( polylog(0.5 + 0.5im, 1.0) )     # the singularity
        @test polylog(1, 2) ≈ -pi*im # comes from -log(complex(1-2))
    
        @testset "     dilogarithm special values for real z" begin
            @test dilog(3.0) == polylog(2, 3.0)

            # for instance, see https://maths.dur.ac.uk/users/herbert.gangl/dilog.pdf (Zagier, 2007)
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

            # https://mathworld.wolfram.com/Dilogarithm.html (from Ramanujan)
            @test polylog(2,  1/3) - (1/6)*polylog(2, 1/9) ≈  pi^2/18 - log(3)^2/6
            @test polylog(2, -1/3) - (1/3)*polylog(2, 1/9) ≈ -pi^2/18 + log(3)^2/6
            @test polylog(2, -1/2) + (1/6)*polylog(2, 1/9) ≈ -pi^2/18 + log(2)*log(3) - log(2)^2/2 - log(3)^2/3
            @test polylog(2,  1/4) + (1/3)*polylog(2, 1/9) ≈  pi^2/18 + 2*log(2)*log(3) - 2*log(2)^2 - 2*log(3)^2/3
            # @test polylog(2, -1/8) +       polylog(2, 1/9) ≈  -log(9/8)^2/2.0 # this one is just outside range
            # Berndt 1994, Gordon and McIntosh 1997, also BBP (2.3)
            @test 36*polylog(2, 1/2) - 36*polylog(2, 1/4) - 12*polylog(2, 1/8) + 6 * polylog(2, 1/64)  ≈  pi^2
            # Lima and Campbell
            @test polylog(2, φ^(-3)) -  polylog(2, -φ^(-3)) ≈  φ^3 * ( pi^2 - 18*log(φ)^2 ) / (3 * (φ^6 - 1) )
            # BBP (2.4)
            @test 4*polylog(2, 1/2) - 6*polylog(2, 1/4) - 2*polylog(2, 1/8) + polylog(2, 1/64)  ≈  log(2)^2
           
            # could add in more of the "polylogarithm ladders" here
            # and there are a few more much more complicated ones in mathworld
            # or general form in BBP (2.16)
         end
        
        @testset "     dilogarithm identities (sometimes called functional equations)" begin
            for i=1:length(Z)
                z = Z[i]
                @test polylog(2, z) + polylog(2, 1/z)   ≈ -pi^2/6.0 - log(Complex(-z))^2/2.0
                @test polylog(2, z) + polylog(2, - z)   ≈ polylog(2, z^2) / 2.0
                @test polylog(2, z) + polylog(2, 1 - z) ≈ pi^2/6.0 - log(Complex(z)) * log(Complex(1-z))

                if real(z) > 0.0
                    # not stated, eg Zagier, but seems to be consistent with Rogers
                    @test polylog(2, -z) - polylog(2, 1 - z) + polylog(2, 1 - z^2)/2.0 ≈ -pi^2/12.0 - log(Complex(z))*log(Complex(z+1))
                end

                # more possibilities at https://functions.wolfram.com/ZetaFunctionsandPolylogarithms/PolyLog/17/ShowAll.html

            
                # Landen’s identity, for z real and <1
                #    Li2(z) + Li2(x/(x-1)) = -(1/2) log(1-x)^2
                 
            end 
        end
        
        @testset "     Spence's function" begin
            @test spence(3.0) == polylog(2, 3.0)
        end
        
        @testset "     Rogers L-function" begin
            # Bytsko, 1999, https://arxiv.org/abs/math-ph/9911012
            #   these are the only five algebraic numbers on [0,1] st output is rational
            @test rogers(0.0)   ≈  0
            @test rogers(0.5)   ≈  1/2 # also from  Rogers (1907), p.189
            @test rogers(1.0)   ≈  1
            @test rogers(ρ)     ≈  3/5
            @test rogers(1-ρ)   ≈  2/5

            @test rogers((3-sqrt(5))/2 ) ≈ 6/15 # from Rogers (1907), p.189

            @test rogers( λ^(-2) ) + rogers( (λ^2-1)^(-2) )      ≈ 4/7 # Bytsko (3.2)
            @test rogers( λ^(-2) ) + rogers( (1+λ)^(-1) )        ≈ 5/7 # Bytsko (3.5)
            @test rogers( sqrt(ρ) ) + rogers( 1/ (1 + sqrt(ρ)) ) ≈ 13/10 # Bytsko (3.16)
            @test rogers( 0.5 - 0.5*ρ) + rogers( 2*ρ -1 )        ≈ 1/2 # Bytsko (3.24)
            @test rogers(1 - 1/sqrt(2)) + rogers(sqrt(2)-1)      ≈ 3/4 # Bytsko (3.14)
            @test rogers(1/sqrt(2)) - rogers(sqrt(2)-1)          ≈ 1/4 # Bytsko (3.14)

            # from Dilogarithm identities and spectra in conformal field theory, Anatol N. Kirillov, arXiv:hep-th/9211137, (0.4)
            @test rogers(ρ^20) - 2*rogers(ρ^10) - 15*rogers(ρ^4) + 10*rogers(ρ^2) ≈ 6/5
            # there are many more in this doc, but none so cute
             
            for i=1:length(Z)
                z = Z[i]
                @test rogers(z) + rogers(1-z)  ≈  1.0 # concise reflection relation

                if real(z) >= 0 # not stated anywhere, but this ID only works for Re(z) >= 0
                    @test rogers(z^2)/2.0  ≈  rogers(z) - rogers( z/(1+z) )   # Abel's duplication formula via Bytsko (3.8),p.6
                    # NB this follows from "Abel's Functional Equation" with L(x) + L(x) (ie x=y)
                    # but that seems to go back to Rogers (1907), (11) and hence (12)
                end
            end
            
            # Zagier "APPENDIX: SPECIAL VALUES AND FUNCTIONAL EQUATIONS OF POLYLOGARITHMS", p.6
            #   but noting that he seems to use the 2nd defintion Gordon and McIntosh (1997) and Loxton (1991, p. 287)
            @test (pi^2/6) * ( 6*rogers(1/3) - rogers(1/9) ) ≈  pi^2 / 3.0

            # there are many more we could put here, but they somewhat double up on the poly checks
        end
        
        @testset "     trilogarithm special values and identities" begin
            @test trilog(3.0) == polylog(3, 3.0)

            # https://mathworld.wolfram.com/Trilogarithm.html
            @test polylog(3,-1.0)             ≈ -3*zeta(3)/4
            @test polylog(3, 0.0)             ≈ 0.0
            @test polylog(3, 0.5)             ≈ log(2)^3/6.0 - pi^2*log(2)/12.0 + (7.0/8.0)*zeta(3)
            @test polylog(3, 1.0)             ≈ zeta(3)
            @test polylog(3, Float64(φ)^(-2)) ≈ 4*zeta(3)/5 + 2*log(φ)^3/3 - 2*pi^2*log(φ)/15 # Duverney, An Introduction to Hypergeometric Functions, p.245, (8.35)

            # https://math.stackexchange.com/questions/932932/known-exact-values-of-the-operatornameli-3-function
            @test polylog(3, 2)            ≈ (pi^2/4) * log(2) + (7/8)*zeta(3) - (pi/2)*log(2)^2*im
            @test polylog(3, Float64(φ)^2) ≈ (4/5)*zeta(3) - (2/3)*log(φ)^3 + (8/15)*pi^2*log(φ) - 2*pi*log(φ)^2*im
            @test polylog(3, im)           ≈ -(3/32)*zeta(3) + (pi^3/32) * im
            @test polylog(3, -im)          ≈ -(3/32)*zeta(3) - (pi^3/32) * im
            @test imag(polylog(3, (1+im)/sqrt(2) )) ≈ 7 * pi^3 / 256
            @test real(polylog(3, (1+im)/2 )) ≈ log(2)^3/48 - (5/192)*pi^2*log(2) + (35/64)*zeta(3)
            @test real(polylog(3, (1+im))) ≈ (pi^2/32)*log(2) + (35/64)*zeta(3)

            for i=1:length(Z)
                z = Z[i]
                # this one isn't accurate enough yet?: @test polylog(3,z) + polylog(3,1-z) + polylog(3,1- 1/z) ≈ zeta(3) + log(z)^2/6 + π^2 * log(z)/6 - log(z)^2*log(1-z)/2

                # Bailey, D. H.; Borwein, P. B.; and Plouffe, S. "On the Rapid Computation of Various Polylogarithmic Constants." Math. Comput. 66, 903-913, 1997.
                @test  36*polylog(3, 1/2) - 18*polylog(3, 1/4) - 4*polylog(3, 1/8) +   polylog(3, 1/64) ≈ (35/2)*zeta(3) - π^2 * log(2)
                @test -24*polylog(3, 1/2) + 18*polylog(3, 1/4) + 4*polylog(3, 1/8) -   polylog(3, 1/64) ≈ 2*log(2)^3 - 7*zeta(3)
                @test (-48*polylog(3, 1/2) + 54*polylog(3, 1/4) +12*polylog(3, 1/8) - 3*polylog(3, 1/64) ≈ 10*log(2)^3 - 2* π^2 *log(2))
            end
       end
        
        @testset "     general s=n case for real z" begin
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
        
        @testset "     tetralogarithm s=4 special values for real z" begin
            @test tetralog(3.0) == polylog(4, 3.0)

            # from https://math.stackexchange.com/questions/1373123/a-conjectured-identity-for-tetralogarithms-operatornameli-4
            @test  96*polylog(4, 1/2) - 54*polylog(4, 1/4) - 8*polylog(4, 1/8) +   polylog(4, 1/64) ≈ 5*log(2)^4 - 2*pi^2*log(2)^2 + 4*pi^4/9

            # there are some others on that list as well. but again, they mostly duplicate polylog tests
            #  but maybe could include more of the ladders
        end
       
        @testset "     general s=n case for complex z" begin
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

    @testset "    particular values |z| == 1 for complex s" begin
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
    
    @testset "    additional identities" begin
        z = 0.5
        for n=1:5
            @test polylog(-n,z) + (-1)^n * polylog(-n, 1/z) ≈ 0.0 
        end

        # for real s, and real z<1, polylog should be real
        S = [-1, 0.1, 2]
        Zr = [-2, -1.0, 0.1, 0.95]
        for i=1:length(S)
            for j=1:length(Zr)
                s = S[i]
                z = Zr[j]
                # println("s = $s; z = $z")
                @test abs( imag( polylog(s,  z) ) ) < 1.0e-12 # actually they usually do better than this
            end
        end

        # for real s, and real z>=1, the imaginary part is given
        S = [-1.5, 0.1, 2]
        Zr = [1.05, 3.0]
        for i=1:length(S)
            for j=1:length(Zr)
                s = S[i]
                z = Zr[j]
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
        S  = [ -1.0 ,  0.5+0.5im,  3.0,   -1.5,         -1.5,        -1.5,      1   ]
        Zc = [ -0.25, 1.0 + im,   -1.0,  -50.0, -400 + 250im, 600 + 600im,  -500.0  ]

        # series 1 cases, polylog(-1, z) ≈ z ./ (1-z).^2
        i = 1
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( z ./ (1-z).^2, 24, 1, 0) )

        # series 2 cases
        i = 2
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 16, 2, 0) )
        
        # series 3 cases,    polylog(3,-1.0)             ≈ -3*zeta(3)/4
        i = 3
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( -3*zeta(3)/4, 36, 3, 0) )
        
        # m-th root multiplication formula
        i = 4
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 36, 22, 0) )
      
        i = 5
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 58, 32, 0) )

        i = 6
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z),  52, 42, 0) )

        # duplication recursion
        i = 7
        s = S[i]
        z = Zc[i]
        @test all( polylog(s, z, Diagnostics()) .≈ ( polylog(s, z), 120, ((3, 3), (3, 3)), 1) )
    end
                   
end

# @testset "Polylogarithm derivatives WRT to s" begin
    # these are a bit made up, as I didn't find to manym but note that
    #    Li_s(1) = ζ(s)     # Riemann zeta function, for Re(s)>1
    # and er know derivatives of Riemann zeta for some cases
    #     https://mathworld.wolfram.com/RiemannZetaFunction.html
    #         z'(2) = (π^2 / 6) * ( γ + log(2*π) - 12 * log(A) ) ~ -0.93754825431...
                                   
    # z = 1.0
    # s = 2.0
    # polylog(s, z) ≈ zeta(s)
    # polylog_ds(s, 0.999999999999999) ≈ (π^2 / 6) * ( γ + log(2*π) - 12 * log(A) ) # NB series doeesn't converge for z=1.0 exactly
    
# end
