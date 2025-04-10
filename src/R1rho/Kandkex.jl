# define a 'safe' square root function that returns 0 for negative inputs
safesqrt(x) = x < 0 ? NaN : sqrt(x)

function estimatekexfromK(K, σΔδ, bf)
    Kinternal = pmean(K) ± pstd(K) # ensure consistent number of particles for calculations here
    dw = 0 ± (σΔδ * 2π * bf)
    kex = safesqrt(Kinternal^2 - dw^2)
    return kex = typeof(kex)(filter(!isnan, kex.particles)) # remove NaNs
end