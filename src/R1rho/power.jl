# Convert power from Hz to dB scale
function convert_Hz_to_dB(νHz, dBref, pref)
    νref = 1 / (4 * pref)
    ratio = νHz / νref

    # Calculate delta_dB
    delta_dB = log10(ratio) * -20

    # Calculate the final power
    return dBref + delta_dB
end
