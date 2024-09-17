function view2d(inputfilename)
    spec = loadnmr(inputfilename)

    dX = data(spec, F1Dim)
    dY = data(spec, F2Dim)
    Xlabel = label(spec, F1Dim)
    Ylabel = label(spec, F2Dim)
    spectitle = label(spec)

    Z = data(spec) / scale(spec)
    σ = spec[:noise] / scale(spec)

    Zdata = Observable(Z)

    ## contour levels
    baselev = exp.(LinRange(0,5,11))
    baselev = sort([baselev; -baselev])
    clev0 = Observable(5σ)
    clev = @lift $clev0 * baselev

    ## prepare figure
    fig = Figure()#resolution=(800,600))

    mainpanel = fig[1,1]
    mainlayout = fig[1, 1] = GridLayout()
    controlpanel = fig[1, 2]
    controllayout = fig[1, 2] = GridLayout()
    # gridcontrols = fig[1, 2]

    # set up contour plot
    contourax = Axis(mainpanel[2,2], xreversed=true, yreversed=true, xlabel=Xlabel, ylabel=Ylabel)

    contourplot = contour!(contourax, dX, dY, Zdata,
        levels=clev,
        colormap=[colorant"red", colorant"black"],
        colorrange=(-0.001σ,0.001σ),
        lowclip=:red,
        highclip=:blue)
        
    xslider = Slider(mainpanel[1,2], range=0:.001:1, startvalue=0.5)
    yslider = Slider(mainpanel[2,1], range=0:.001:1, startvalue=0.5, horizontal=false)

    # set up control panel
    □contourup = Button(fig, label="contours ↑")
    □contourdown = Button(fig, label="contours ↓")

    □xslice = Toggle(fig, active=false)
    □yslice = Toggle(fig, active=false)
    # labels = [Label(fig, "x slice"), Label(fig, "y slice")]


    # fig[1, 2] = grid!(hcat(toggles, labels), tellheight = false)


    controllayout[1,1] = □contourup
    controllayout[2,1] = □contourdown
    controllayout[3,1] = grid!([Label(fig, "x slice") ;; □xslice])
    controllayout[4,1] = grid!([Label(fig, "y slice") ;; □yslice])

    controllayout.tellheight[] = false


    # contour logic
    on(□contourup.clicks) do clicks
        clev0[] *= 1.3
    end

    on(□contourdown.clicks) do clicks
        clev0[] /= 1.3
    end

    # slice logic
    xsliderpos = lift(contourax.finallimits, xslider.value) do lims, ix
        x0 = lims.origin[1]
        dx = lims.widths[1]
        x = x0 + (1-ix)*dx  # (1-ix) because axis is reversed
        nearest(dX, x)
    end
    ysliderpos = lift(contourax.finallimits, yslider.value) do lims, iy
        y0 = lims.origin[2]
        dy = lims.widths[2]
        y = y0 + (1-iy)*dy
        nearest(dY, y)
    end
    yslicedata = lift(contourax.finallimits, xsliderpos, Zdata) do lims, xpos, Z
        x0 = lims.origin[1]
        x1 = x0 + lims.widths[1]
        xmid = (x0 + x1) / 2
        dx = lims.widths[1] / 4

        y0 = lims.origin[2]
        y1 = y0 + lims.widths[2]

        ix = findnearest(dX, xpos)
        iy0 = findnearest(dY, y0)
        iy1 = findnearest(dY, y1)
        iy = min(iy0, iy1):max(iy0, iy1)

        slice = Z[ix, iy]
        zscale = dx / maximum(abs.(slice))
        xs = xmid .- slice*zscale
        ys = dY[iy]

        Point2f[zip(xs, ys)...]
    end
    xslicedata = lift(contourax.finallimits, ysliderpos, Zdata) do lims, ypos, Z
        y0 = lims.origin[2]
        y1 = y0 + lims.widths[2]
        ymid = (y0 + y1) / 2
        dy = lims.widths[2] / 4

        x0 = lims.origin[1]
        x1 = x0 + lims.widths[1]

        iy = findnearest(dY, ypos)
        ix0 = findnearest(dX, x0)
        ix1 = findnearest(dX, x1)
        ix = min(ix0, ix1):max(ix0, ix1)

        slice = vec(Z[ix, iy])
        zscale = dy / maximum(abs.(slice))
        ys = ymid .- slice*zscale
        xs = dX[ix]

        Point2f[zip(xs, ys)...]
    end

    xline = hlines!(contourax, ysliderpos, color=:deepskyblue, linestyle=:dash)
    xslice = lines!(contourax, xslicedata, color=:deepskyblue)
    yline = vlines!(contourax, xsliderpos, color=:chartreuse3, linestyle=:dash)
    yslice = lines!(contourax, yslicedata, color=:chartreuse3)

    connect!(xline.visible, □xslice.active)
    connect!(yline.visible, □yslice.active)
    connect!(xslice.visible, □xslice.active)
    connect!(yslice.visible, □yslice.active)

    # display figure
    fig
end