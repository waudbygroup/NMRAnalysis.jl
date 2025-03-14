function setupR1rhopowers()
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
    println("$magenta Enter p1 (in us):$reset")
    print("> ")
    p1 = parse(Float64, readline()) * 1e-6 # convert to seconds
    
    println()
    println("$magenta Input your own list of target spinlock strengths (in Hz) separated by commas, or press enter for a default list:$reset")
    print("> ")
    input = readline()
    
    if input == ""
        target_spinlock_strengths = [100, 200, 300, 500, 750, 1000, 1500, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000]
    else
        target_spinlock_strengths = parse.(Float64, split(input, ","))
    end
    
    # Calculate the final powers
    final_powers = convert_Hz_to_dB.(target_spinlock_strengths, pldb1, p1)
    
    # Shuffle the final powers list and the corresponding spinlock strengths
    shuffled_indices = shuffle(1:length(final_powers))
    shuffled_final_powers = final_powers[shuffled_indices]
    shuffled_spinlock_strengths = target_spinlock_strengths[shuffled_indices]
    
    # Print the final powers in the specified format
    println()
    println("The default spinlock strengths are:", shuffled_spinlock_strengths)
    println()
    println("Copy & paste the list provided between the dashed lines.")
    println(line_break)
    println("dB")
    for power in shuffled_final_powers
        println(@sprintf("%.2f", power))
    end
    println(line_break)
end