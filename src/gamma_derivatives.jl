using SpecialFunctions
using Base.MathConstants
using DataFrames, CSV
using Printf

"""
    gamma_derivatives()

Derive a table of derivatives of the Gamma function at 1, i.e., Gamma^{(m)}(a)

This code is used to create the results, which are included in the main code via a lookup table.

"""
function gamma_derivatives(k::Integer, x::Float64)
    #    general form from: http://erikerlandson.github.io/blog/2016/06/15/computing-derivatives-of-the-gamma-function/
    return SpecialFunctions.gamma(x) * D(k,0,x)
end

function D(k::Integer, n::Integer, x::Float64)
    # general form from: http://erikerlandson.github.io/blog/2016/06/15/computing-derivatives-of-the-gamma-function/
    if n<0
        throw(DomainError(n))
    elseif k<0
        throw(DomainError(k))
    end
    if k==0
        return 1
    elseif k==1
        return SpecialFunctions.polygamma(n,x)
    else
        total = D(k-1,n+1,x)
        for j=0:n
            total += binomial(n,j) * D(1,j,x) * D(k-1,n-j,x)
        end
        return total
    end
end

function g1(t::Int) # Crandall,2012, p.17
    # t derivate of Gamma function at 1
    PG21 = -2.4041138063191902 # SpecialFunctions.polygamma(2,1)
    PG41 = -24.8862661234409   # SpecialFunctions.polygamma(4,1)
    if t==0
        return 1.0
    elseif t==1
        return -γ
    elseif t==2
        # return γ^2 - γ + pi^2/6, Crandall seems to be wrong on this
        return γ^2 + pi^2/6  # from Mathematica
    elseif t==3
        # return -2*γ^3 + 9*γ^2 - (π^2+6)*γ + 3*π^2/2 - 4*SpecialFunctions.zeta(3), Crandall seems to be wrong on this
        return -γ^3  - π^2*γ/2 + PG21 # from Mathematica
    elseif t==4
        return γ^4 + γ^2*π^2 + 3*π^4/20 - 4*γ*PG21 # from Mathematica
    elseif t==5
        return -γ^5 - (20/12)*γ^3*π^2 - (9/12)*γ*π^4 + (10*γ^2 + (20/12)*π^2)*PG21 + PG41 # from Mathematica
    else
        return NaN
    end
end

function g2(t::Int) # http://erikerlandson.github.io/blog/2016/06/15/computing-derivatives-of-the-gamma-function/
    # t derivate of Gamma function at 1
    return gamma_derivatives(t,1.0) 
end


# #### output a table
# K = 14
# k = collect(0:K)
# G1 = zeros(Float64, K+1)
# G2 = zeros(Float64, K+1)
# for i=1:K+1
#     G1[i] = g1(i-1)
#     G2[i] = g2(i-1)
#     @printf("G1[i] = %.20f, G2[i] = %.20f\n", G1[i], G2[i])
# end

# df = DataFrame( k=k, G1=G1, G2=G2 )
# CSV.write("../data/gamma_derivatives.csv", df)

