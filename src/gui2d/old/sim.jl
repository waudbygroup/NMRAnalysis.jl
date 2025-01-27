function simZ!(state)
    specdata = state[:specdata]
    z = specdata[:zfit]
    for slice in 1:specdata[:nspec]
        z[slice] *= 0
        xaxis = dims(specdata[:spectra][slice], F1Dim)
        yaxis = dims(specdata[:spectra][slice], F2Dim)
        for peak in state[:peaks][]
            position = peak.sim_parameters[:position][slice]
            isnothing(position) && continue
            amplitude = peak.sim_parameters[:amplitude][slice]
            R2X = peak.sim_parameters[:R2X][slice]
            R2Y = peak.sim_parameters[:R2Y][slice]
            simx = lineshape(xaxis, position[1], R2X)
            simy = lineshape(yaxis, position[2], R2Y)
            z[slice] += amplitude * simx .* simy'
        end
    end
end