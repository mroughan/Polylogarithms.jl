include("test_defs.jl")
include("../src/zeta_derivative.jl")

@testset "Riemann zeta derivatives" begin

#    @testset "    throws errors" begin
#        @test_throws DomainError stieltjes(-1)
#        @test_throws DomainError stieltjes(79)
#        @test_throws MethodError stieltjes(1.5)
#    end
    
    @testset "    types" begin
        @test typeof( zeta_derivative(2.0)[1] ) == Float64
        @test typeof( zeta_derivative(2.0 + 1.0*im)[1] ) == Complex{Float64}
    end

    @testset "    values" begin
        # https://mathworld.wolfram.com/RiemannZetaFunction.html
        # https://en.wikipedia.org/wiki/Particular_values_of_the_Riemann_zeta_function#Derivatives
        # Also COMPUTATION OF THE DERIVATIVES OF THE RIEMANN ZETA-FUNCTION IN THE COMPLEX DOMAIN, U.S.S.R. Comput.Maths.Math.Phys.,Vo1.28,No.4,pp.115-124,1988, A.YU. YEREMIN, I.E. KAPORIN and M.K. KERIMOV
        #       but the value for (1/2) seems different
        
        @test zeta_derivative(2)[1] ≅ (1/6) * π^2 * ( γ + log(2π) - 12*log(A) )
                                  # = -0.93754825431...
        
        @test zeta_derivative(1/2)[1] ≅ ( π + 2*γ + 6*log(2) + 2*log(π) ) * zeta(1/2) / 4
                # = -3.92264613920915172747153144671459951373032397150650...
                # also see https://oeis.org/A114875

        @test zeta_derivative(-2.0)[1] ≅ test_d_near_zeros( -2 )
        @test zeta_derivative(-4.0)[1] ≅ test_d_near_zeros( -4 )
        @test zeta_derivative(-6.0)[1] ≅ test_d_near_zeros( -6 )
        @test zeta_derivative(-8.0)[1] ≅ test_d_near_zeros( -8 )
        # could use general form for ζ(-2n), see https://mathworld.wolfram.com/RiemannZetaFunction.html
        
        @test zeta_derivative(-1)[1] ≅ (1/12) - log(A)

        
        # could do expantion around zero better
        #  very simple to do as a special case, but should do expansion around zero
        @test zeta_derivative(0.0)[1] ≈ -(1/2) * log(2π) # not quite meeting high accuracy goal 

        # from https://en.wikipedia.org/wiki/Particular_values_of_the_Riemann_zeta_function#Derivatives
        @test zeta_derivative(3)[1] ≅ −0.19812624288563685333
        @test zeta_derivative(-1/2)[1] ≅ −0.36085433959994760734
        @test zeta_derivative(-3)[1] ≅ +0.005378576357774301144
        @test zeta_derivative(-5)[1] ≅ −0.00057298598019863520499
        @test zeta_derivative(-7)[1] ≅ −0.00072864268015924065246
        
    end

end

stop

# look at a plot near zero
S= 10 .^ collect( -15 : 0.1 : 0.5 )

ze = [ zeta_derivative.(s)[1] for s in S]
ze2 = [ zeta_series_4.(s)[1] for s in S]
yb = 1.0e-2
p1 = plot( log10.(S), ze .+ log(2π)/2 ; ylim=[-yb, yb] )
plot!( log10.(S), ze2 .+ log(2π)/2 )

ze3 = [ zeta_derivative.(-s)[1] for s in S]
ze4 = [ zeta_series_4.(-s)[1] for s in S]
plot!(p1, log10.(S), ze3 .+ log(2π)/2 )
plot!(p1, log10.(S), ze4 .+ log(2π)/2 )





# near -2
m = -6
S= m .+ 10 .^ collect( -17 : 0.1 : -0.1 )

ze = [ zeta_reflection.(s)[1] for s in S]
ze2 = [ zeta_series_4.(s)[1] for s in S]
ze3 = [ -zeta(s)/(s-m) for s in S]
ze4 = [ 2^(s-1) * pi^s * cos(pi*s/2) * gamma(1-s) * zeta(1-s) for s in S]
ze5 = [ zeta_d_near_zeros.(s)[1] for s in S]
yb = 1.0e-4
p2 = plot( log10.(S .- m), ze .- test_d_near_zeros(m) ; ylim=[-yb, yb] )
plot!(p2, log10.(S .- m), ze2 .- test_d_near_zeros(m) )
plot!(p2, log10.(S .- m), ze3 .- test_d_near_zeros(m) )
plot!(p2, log10.(S .- m), ze4 .- test_d_near_zeros(m) )
plot!(p2, log10.(S .- m), ze5 .- test_d_near_zeros(m) )


p3 = plot(  log10.(S .- m), log10.(abs.( ze .- test_d_near_zeros(m) )) ; label="reflection")
plot!(p3,   log10.(S .- m), log10.(abs.( ze2 .- test_d_near_zeros(m) )) ; label="alternating series" )
plot!(p3,   log10.(S .- m), log10.(abs.( ze4 .- test_d_near_zeros(m) )) ; label="main term" )
plot!(p3,   log10.(S .- m), log10.(abs.( ze5 .- test_d_near_zeros(m) )) ; label="final form")
plot!( xlabel="log10(s - m)", ylabel="log10(abs(error))" )

# tables of values from
#   WALTHER A., Anschauliches zur Riemannschen Zeta-Function.Acta Math., 46, 393-400, 1926.
#   Via COMPUTATION OF THE DERIVATIVES OF THE RIEMANN ZETA-FUNCTION IN THE COMPLEX DOMAIN, U.S.S.R. Comput.Maths.Math.Phys.,Vo1.28,No.4,pp.115-124,1988, A.YU. YEREMIN, I.E. KAPORIN and M.K. KERIMOV
#


# comparison to Mathematica or mpmath
#   
