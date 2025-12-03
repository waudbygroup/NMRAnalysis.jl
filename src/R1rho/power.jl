# Convert power from Hz to dB scale
function convert_Hz_to_dB(νHz, powerref, p90ref)
    νref = 1 / (4 * p90ref)
    ratio = νHz / νref

    # Calculate delta_dB
    delta_dB = log10(ratio) * -20

    # Calculate the final power
    return db(powerref) + delta_dB
end
