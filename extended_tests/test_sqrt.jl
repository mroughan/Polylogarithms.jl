using StatsBase
using PyPlot

function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

x = -5 : 0.13 : 5
y = -5 : 0.13 : 5
(X,Y) = meshgrid(x, y)
Z = X + Y*im

Mu1 = log.(Z)./(2*pi)
Mu2 = log.(-sqrt.(Z))./(2*pi)
Mu3 = log.(sqrt.(sqrt.(Z)))./(2*pi)

mean( abs.(Mu2) ./ abs.(Mu1))


theta = -pi : 0.01 : pi
r = 0.5
z = r .* exp.(im*theta)
z2 = sqrt.(z)
z3 = -z2

figure()
plot(theta,  real.(z), imag.(z) )
plot(theta,  real.(z2), imag.(z2), "--" )
axis("equal")

figure()
plot(theta, abs.(log.(z)) )
plot(theta,  abs.(log.(z2)), "--"  )
plot(theta,  abs.(log.(z3)), ":"  )

