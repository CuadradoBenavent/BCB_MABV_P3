#!/bin/bash
#Input file name
input_file="trajectory-6.xyz"

# Output file name
output_file="energies-6.dat"

# Use awk to extract iteration and energies
awk 'BEGIN{print "Iteration", "Potential Energy", "Kinetic Energy", "Total Energy"} /Iteration:/ {iteration = $2} /Potential/ {print iteration, $3, $6, $9}
    ' "$input_file" > "$output_file"

    echo "Energies extracted and saved to $output_file"
