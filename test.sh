#!/bin/bash

# Set the depth level (n)
depth=$1
maxPatchSize=$2
execPath=$3

# Define the base paths for models
model_base_path="meshes/mambo"

# Initialize variables for accumulating results
total=0
count=0

# Function to run the test and check the result
run_test() {
    model_path="$1"

    echo "Running test with model: $model_path"
    
    # Run the program and capture its output
    $execPath model=$model_path result_path=output/ maxPatchSize=$maxPatchSize
    result=$?

    # Accumulate the result
    total=$(echo "$total + $result" | bc)
    count=$((count + 1))

    if (( $(echo "$result > 100" | bc -l)  )); then
        echo "Test failed with result: $result"
        exit 1
    fi

    # Check if the result is superior to 60
    if (( $(echo "$result <= 15" | bc -l) )); then

        echo "Test failed with result: $result"
        exit 1
    else
        echo "Test passed with result: $result"
    fi
}

# Loop through the depth levels and run the tests
for ((i=1; i<=depth; i++)); do
    for model in Basic/B Simple/S Medium/M; do
        model_file="$model_base_path/$model${i}.mesh"
        run_test "$model_file"
    done
done

# Calculate the average result
average=$(echo "scale=2; $total / $count" | bc)
echo "All tests passed successfully. Average result: $average"
