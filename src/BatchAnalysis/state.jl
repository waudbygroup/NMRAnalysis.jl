function initialisestate(spectra; referenceshift=0.0, standardx=0.0, standarddx=0.05, unknownx=7.0, unknowndx=0.05)
    state = Dict{String, Any}()

    # spectra
    state["spectra"] = Observable(spectra)
    state["nspectra"] = lift(length, state["spectra"])
    state["labels"] = lift(label, state["spectra"])
    state["selected"] = Observable(fill(true, state["nspectra"][]))
    state["nselected"] = lift(sum, state["selected"])
    state["selectedspectra"] = lift((spectra, selected) -> spectra[selected], state["spectra"], state["selected"])

    # chemical shift referencing
    state["referenced"] = Observable(false)
    state["referenceshift"] = Observable(referenceshift)

    # integration
    state["integrating"] = Observable(false)
    state["standardx"] = Observable(standardx)
    state["standarddx"] = Observable(standarddx)
    state["unknownx"] = Observable(unknownx)
    state["unknowndx"] = Observable(unknowndx)

    state["standardx1"] = lift((x,dx) -> x - dx/2, state["standardx"], state["standarddx"])
    state["standardx2"] = lift((x,dx) -> x + dx/2, state["standardx"], state["standarddx"])
    state["unknownx1"] = lift((x,dx) -> x - dx/2, state["unknownx"], state["unknowndx"])
    state["unknownx2"] = lift((x,dx) -> x + dx/2, state["unknownx"], state["unknowndx"])

    # filenames
    state["figurefilename"] = Observable("spectra.pdf")
    state["resultsfilename"] = Observable("results.csv")

    return state
end


function preparegui!(state)
    state["gui"] = Dict{String, Any}()
    gui = state["gui"]

    gui["plotscale"] = Observable(1.0)
    gui["linedata"] = lift(linedata, state["selectedspectra"], gui["plotscale"])
    gui["linecolors"] = lift(linecolors, state["nselected"])
end