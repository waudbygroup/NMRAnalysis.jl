function analyseAK(A, K)
    pBs = logrange(.0001, .02, 100)

    # Initialize empty arrays to store all k and pB pairs
    all_ks = Float64[]
    all_pBs = Float64[]
    all_dws = Float64[]

    for pB in pBs
        # Call your solver function
        ks = solveAKfork(A, K, pB)
        dws = map(k->sqrt(K^2-k^2), ks)
        
        # For each k value returned, add a (pB, k) pair
        for (k,dw) in zip(ks,dws)
            push!(all_ks, k)
            push!(all_pBs, pB)
            push!(all_dws, dw)
        end
    end

    f = Figure()
    ax = Axis(f[1,1], ylabel="koff and Δω (s-1)", xlabel="p_bound")
    scatter!(ax, all_pBs, all_ks, label="koff")
    scatter!(ax, all_pBs, all_dws, label="Δω")
    axislegend(ax;position=:lt)
    xlims!(0, 0.0205)
    ylims!(0, K)
    
    ax2 = Axis3(f[1,2], xlabel="p_bound", ylabel="koff", zlabel="Δω", aspect=:equal)
    scatter!(ax2,all_pBs, all_ks, all_dws)

    f
end

function solveAKfork(A, K, pB=0.01)
    # k^3 - K^2*k + A*K^2/pB = 0
    p = -K^2
    q = A*K^2/pB
    t1 = (-q + sqrt(Complex(q^2 + 4(p/3)^3))) / 2
    t2 = (-q - sqrt(Complex(q^2 + 4(p/3)^3))) / 2
    u = t1^(1/3)
    v = t2^(1/3)
    x1 = u + v
    ω = -0.5 + sqrt(3)/2im
    x2 = ω*u + ω^2*v
    x3 = ω^2*u + ω*v
    return clean_and_filter_solutions([x1, x2, x3])
end

function clean_and_filter_solutions(solutions; tol=1e-10)
    # Initialize empty array for the filtered results
    filtered_solutions = Float64[]
    
    for z in solutions
        # If imaginary part is relatively small compared to real part or absolutely tiny
        if abs(imag(z)) < tol * max(1.0, abs(real(z))) || abs(imag(z)) < tol
            # Convert to real
            real_value = real(z)
            
            # Only keep positive values
            if real_value > 0
                push!(filtered_solutions, real_value)
            end
        end
    end
    
    return filtered_solutions
end