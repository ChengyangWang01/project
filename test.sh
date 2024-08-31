#!/bin/bash

# Program names
program1="main"
program2="cylinder"

# Output directory
output_dir="outcome"

# Ensure the output directory exists
mkdir -p $output_dir

# List of process counts for which the programs will be run
process_counts=(1 2 4 8 16 32 64)

# Loop over each process count
for np in "${process_counts[@]}"; do
    # Running first program
    echo "Running $program1 with $np processes:"
    (cd build && mpirun -np $np ./$program1 > "../$output_dir/np${np}_${program1}.txt") || {
        echo "Failed to run $program1 with $np processes."
        continue
    }

    # Running second program
    echo "Running $program2 with $np processes:"
    (cd build && mpirun -np $np ./$program2 > "../$output_dir/np${np}_${program2}.txt") || {
        echo "Failed to run $program2 with $np processes."
        continue
    }
done
