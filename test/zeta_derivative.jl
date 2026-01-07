include("test_defs.jl")

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
        # and https://oeis.org/search?q=derivative+of+Riemann%27s+zeta&go=Search
        
        @test zeta_derivative(2)[1] ≅ (1/6) * π^2 * ( γ + log(2π) - 12*log(A) )
                                  # = -0.93754825431...
        
        @test zeta_derivative(1/2)[1] ≅ ( π + 2*γ + 6*log(2) + 2*log(π) ) * zeta(1/2) / 4
                # = -3.92264613920915172747153144671459951373032397150650...
                # also see https://oeis.org/A114875

        @test zeta_derivative(-2.0)[1] ≅ test_d_near_zeros( -2 )
        @test zeta_derivative(-4.0)[1] ≅ test_d_near_zeros( -4 )
        @test zeta_derivative(-6.0)[1] ≅ test_d_near_zeros( -6 )
        @test zeta_derivative(-8.0)[1] ≅ test_d_near_zeros( -8 )
        @test zeta_derivative(-10.0)[1] ≅ test_d_near_zeros( -10 )
        @test zeta_derivative(-12.0)[1] ≅ test_d_near_zeros( -12 )
        # could use general form for ζ(-2n), see https://mathworld.wolfram.com/RiemannZetaFunction.html
        
        @test zeta_derivative(-1)[1] ≅ (1/12) - log(A)

        # could do expantion around zero better
        #  very simple to do as a special case, but should do expansion around zero
        @test zeta_derivative(0.0)[1] ≈ -(1/2) * log(2π) # not quite meeting high accuracy goal 

        # from https://en.wikipedia.org/wiki/Particular_values_of_the_Riemann_zeta_function#Derivatives
        @test zeta_derivative(3)[1] ≅ −0.19812624288563685333     # OEIS: A244115
        @test zeta_derivative(-1/2)[1] ≅ −0.36085433959994760734  # OEIS: A271854
        @test zeta_derivative(-3)[1] ≅ +0.005378576357774301144   # OEIS: A259068
        @test zeta_derivative(-5)[1] ≅ −0.00057298598019863520499 # OEIS: A259070


        # we are loosing a bit of accuracy arounf larger negative values, so use less strict criteria
        @test zeta_derivative(-7)[1] ≈ -0.000728642680159240652467233354650360611902887720925418318636386154 # OEIS: A259072, https://oeis.org/A259072
        # zeta_derivative(-7)[1] ≅ -121/11200 + (γ + log(2π))/240 - 315*zeta_derivative(8)[1]/(8π^8)

        @test zeta_derivative(-9)[1] ≈ 0.0031301453197885727549257682907854467026693658654815 # A266260
        @test zeta_derivative(-11)[1] ≈ -0.012752984479966656113522525488725798156238937049874292793246366661 # A266262
        @test zeta_derivative(-13)[1] ≈  0.06374987374457688028603868147333505564882735 # A260660
        @test zeta_derivative(-15)[1] ≈ -0.400319302807725593843580317520320367201261286266232944284106942 # A266270
        @test zeta_derivative(-17)[1] ≅ test_d_near_zeros(-17)

        # I can't find any identities that apply for complex arguments??? 


        # zeros from "The Riemann Zeta-Function and Its Derivatives", Bejoy K. Choudhury, Choudhury-RiemannZetaFunctionDerivatives-1995.pdf
        #   Table 1 (but note that the zeros aren't exact, and there are 10 listed, so I can add
        zeta_derivative_zeros = [
#            -2.712628292045741016 
            -4.9367621085949478689 
        ]
        for i=1:size(zeta_derivative_zeros, 1)
            @test zeta_derivative( zeta_derivative_zeros[i,1] )[1]  ≅ 0.0
        end
    end
end

