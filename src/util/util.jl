module Util

include("maybevector.jl")

export MaybeVector, SingleElementVector, StandardVector
export nearest, findnearest
export extract_residue_number


nearest(A::AbstractArray, t) = findmin(abs.(A .- t))[1]
findnearest(A::AbstractArray, t) = findmin(abs.(A .- t))[2]

"""
    extract_residue_number(label::AbstractString)::Int

Convert a residue label to its corresponding residue number. Handles various formats:
- Standard residue labels (e.g., "A7", "D98", "G321")
- Reversed format (e.g., "7A", "98D", "321G")
- Methyl groups (e.g., "I13CD1", "L98CD2", "M98CE")
- Non-standard residues (e.g., "X99" -> -99)

Returns:
- The residue number as an integer
- For non-standard residues (not in standard_residues), returns negative number
"""
function extract_residue_number(label)
    # Define standard amino acid one-letter codes
    standard_residues = Set(['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                            'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'])
    
    # Handle empty input
    isempty(label) && throw(ArgumentError("Empty label provided"))
    
    # Extract all digits from the label
    numbers = match(r"\d+", label)
    
    # Return 0 if no numbers found
    if isnothing(numbers)
        return 0
    end
    
    residue_num = parse(Int, numbers.match)
    
    # Find the first letter in the label
    first_letter = first(filter(isletter, label))
    
    # Check if it's a standard residue
    if !in(first_letter, standard_residues)
        residue_num = -residue_num
    end
    
    return residue_num
end

end