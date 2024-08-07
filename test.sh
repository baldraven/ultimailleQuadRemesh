#!/bin/bash

# Set the depth level (n)
depth=$1
execPath=$2
extraParams=$3

model_base_path="meshes/mambo"

total=0
count=0

start_time=$(date +%s)

run_test() {
    model_path="$1"

    echo "Running test with model: $model_path"
    
    # Run the program and capture its output
    #$execPath model=$model_path $extraParams
    $execPath model=$model_path result_path=output/ $extraParams 

    result=$?

    total=$(echo "$total + $result" | bc)
    count=$((count + 1))

    if (( $(echo "$result > 100" | bc -l)  )); then
        echo "Test failed with result: $result"
        exit 1
    fi

    if (( $(echo "$result <= 5" | bc -l) )); then

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

end_time=$(date +%s)
# Calculate the average result
average=$(echo "scale=2; $total / $count" | bc)

# Calculate the time spent
time_spent=$((end_time - start_time))
average_time_per_3_meshes=$(echo "scale=2; $time_spent / $depth" | bc)

echo "\nAll tests passed successfully. Average result: $average"
echo "Time spent: $time_spent seconds (Average per 3 meshes: $average_time_per_3_meshes seconds)"
