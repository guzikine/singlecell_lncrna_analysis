#!/bin/bash

# Check if input file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <input.fastq>"
    exit 1
fi

# Input FASTQ file
INPUT_FASTQ=$1

# Calculate the average read length
AVG_READ_LENGTH=$(gzip -dc $INPUT_FASTQ | awk '{if(NR%4==2) print length($0)}' | awk '{s+=$1} END {print s/NR}')

# Calculate the standard deviation of read lengths
READ_SD=$(gzip -dc $INPUT_FASTQ | awk '{if(NR%4==2) lengths[NR/4] = length($0)} END { 
    # Calculate mean
    for (i in lengths) sum += lengths[i];
    mean = sum / length(lengths);
    
    # Calculate variance
    for (i in lengths) sum_sq += (lengths[i] - mean) ^ 2;
    
    # Calculate standard deviation
    print sqrt(sum_sq / length(lengths))
}')

READ_SD=$(echo "$READ_SD + 0.1" | bc)

# Output the results
echo "$AVG_READ_LENGTH $READ_SD"