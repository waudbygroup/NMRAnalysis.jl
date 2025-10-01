"""
    convert_valist_to_Hz(expt)

Converts a valist to spinlock strengths in Hz.

# Arguments
- `expt`: The experiment with a valist to be converted.

# Returns
An array of spinlock strengths, in Hz.
"""
function convert_valist_to_Hz(expt)
    # Extract valist from expt
    valist = acqus(expt, :valist)
    # Convert to Hz
    return hz.(valist, acqus(expt, :pl, 1), acqus(expt, :p, 1), 90)
end

function convert_Hz_to_dB(νHz, dBref, pref)
    νref = 1 / (4 * pref)
    ratio = νHz / νref

    # Calculate delta_dB
    delta_dB = log10(ratio) * -20

    # Calculate the final power
    return dBref + delta_dB
end
