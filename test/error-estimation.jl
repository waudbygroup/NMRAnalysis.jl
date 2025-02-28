using GLMakie

koff = LinRange(100, 10000, 200)
dw = LinRange(100, 3500*1.5, 200)
dd = dw ./ (2π*565) # 600 MHz 1H

safesqrt = x -> x < 0 ? NaN : sqrt(x)
relerr(koff, dw, dwav) = 100*(safesqrt(1+(dw^2-dwav^2)/(koff^2))-1)

z = [relerr(koff, dw, 3500) for koff in koff, dw in dw]

f = Figure()
ax = Axis(f[1, 1], xlabel="koff / s-1", ylabel="Δω / s-1")
hm = heatmap!(ax, koff, dw, z, colorrange=(-50,50), colormap=:bluesreds)

Colorbar(f[1,2], hm, label="relative error (%)")

f

##
ddav = 1
dwav = ddav * 2π * 565 # 600 MHz 1H
f = Figure()
ax = Axis(f[1, 1], xlabel="True Δδ / ppm", ylabel="Error in dissociation rate (%)")
vlines!(ax, [ddav], color=:black, label="Assumed Δδ")
lines!(ax, dd, [abs.(relerr.(2500, dw, dwav)) for dw in dw], label="B = 2,500 s-1")
lines!(ax, dd, [abs.(relerr.(5000, dw, dwav)) for dw in dw], label="B = 5,000 s-1")
lines!(ax, dd, [abs.(relerr.(10000, dw, dwav)) for dw in dw], label="B = 10,000 s-1")
lines!(ax, dd, [abs.(relerr.(20000, dw, dwav)) for dw in dw], label="B = 20,000 s-1")
axislegend(;position=:lt)
f