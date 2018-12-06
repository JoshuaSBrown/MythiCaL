#!/bin/bash

sigma=0.07
length=11
particles=3
threshold=1000
time=0.001
total_crude=0
total_coarse=0
echo "Crude    Coarse" > time_all.txt
number_of_seeds=10
for seed in $(seq 1 $number_of_seeds)
do
  ./performance_test_crude_vs_coarsegrain $sigma $length $seed  $particles $threshold $time > out.txt
  value_crude=$(cat out.txt | grep "Crude Monte Carlo Run Time" | awk '{print $6}') 
  total_crude=$(( $value_crude+$total_crude ))
  value_coarse=$(cat out.txt | grep "Coarse Monte Carlo Run Time" | awk '{print $6}') 
  total_coarse=$(( $value_coarse+$total_coarse )) 
  echo $value_crude"          "$value_coarse >> time_all.txt
done

average_time_crude=$(bc <<< "scale = 6; ( $total_crude / $number_of_seeds )")
average_time_coarse=$(bc <<< "scale = 6; ( $total_coarse / $number_of_seeds )")

echo $average_time_crude
echo $average_time_coarse

echo $sigma" "$length" "$particles" "$average_time_crude >> average_crude_time.txt
echo $sigma" "$length" "$particles" "$average_time_coarse >> average_coarse_time.txt
