"""
    convert_valist_to_Hz(expt)

Converts a valist to spinlock strengths in Hz.

# Arguments
- `expt`: The experiment with a valist to be converted.

# Returns
An array of spinlock strengths, in Hz.
"""
function convert_valist_to_Hz(expt)
    # Extract valist from expt (in dB)
    valist = acqus(expt, :valist)
    
    return convert_dB_to_Hz(valist, expt)
end

function convert_dB_to_Hz(dB, expt)
    # Get the reference power level
    plw1 = acqus(expt, :plw, 1) # pl1 in W
    dBref = convert_W_to_dB(plw1) # convert pl1 from W to dB
    
    # Get the reference pulse length
    p1 = acqus(expt, :p, 1) * 1e-6 # p1 in seconds 
    νref = 1/(4*p1) # p1 in seconds 

    # Apply the conversion to the input power
    ΔdB = dB .- dBref
    return @. νref * 10^(-ΔdB / 20)
end

convert_W_to_Hz(plW, expt) = convert_dB_to_Hz(convert_W_to_dB(plW), expt)

function convert_W_to_dB(plW)
    return -10 * log10(plW)
end