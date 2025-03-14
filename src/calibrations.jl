using Printf
using Random

# Function to calculate the final powers
function calculate_final_powers(p1, pldb1, target_spinlock_strengths)
    # Calculate the power in Hz
    power = 1 / (p1 * 4) * 10^6
    
    # Initialize an empty list to store the final powers
    final_powers = []
    
    # Iterate over each target spinlock strength
    for target_spinlock in target_spinlock_strengths
        # Calculate the ratio
        ratio = target_spinlock / power
        
        # Calculate delta_dB
        delta_dB = log10(ratio) * -20
        
        # Calculate the final power
        final_power = pldb1 + delta_dB
        
        # Round the final power to the hundredths place
        rounded_final_power = round(final_power, digits=2)
        
        # Append the rounded final power to the list
        push!(final_powers, rounded_final_power)
    end
    
    return final_powers
end

# Main script
function main()
    # ANSI escape code for magenta
    magenta = "\033[35m"
    reset = "\033[0m"
    
    # Get the terminal width
    term_width = displaysize(stdout)[2]
    line_break = repeat("-", term_width)
    
    # Prompt the user for inputs with magenta color and caret
    println()
    println("$magenta Enter pldb1:$reset")
    print("> ")
    pldb1 = parse(Float64, readline())
    
    println()
    println("$magenta Enter p1:$reset")
    print("> ")
    p1 = parse(Float64, readline())
    
    println()
    println("$magenta Input your own list of target spinlock strengths (in Hz) separated by commas, or type 'd' for a default list:$reset")
    print("> ")
    input = readline()
    
    if lowercase(input) == "d"
        target_spinlock_strengths = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    else
        target_spinlock_strengths = parse.(Float64, split(input, ","))
    end
    
    # Calculate the final powers
    final_powers = calculate_final_powers(p1, pldb1, target_spinlock_strengths)
    
    # Shuffle the final powers list
    shuffled_final_powers = shuffle(final_powers)
    
    # Print the final powers in the specified format
    println()
    println("Copy & paste the list provided between the dashed lines.")
    println(line_break)
    println("dB")
    for power in shuffled_final_powers
        println(@sprintf("%.2f", power))
    end
    println(line_break)
end

# Run the main script
main()