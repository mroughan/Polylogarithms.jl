using PyPlot

s = im
z = 0.5
n = 20
k = 1:n
w = z.^k ./ k.^s

figure(1)
plot( real(w), imag(w), "o-")
axis("equal")
xlim([-0.5, 0.5])
