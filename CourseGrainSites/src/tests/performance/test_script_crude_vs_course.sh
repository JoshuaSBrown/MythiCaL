#!/bin/bash

resolution=10
sigma=0.07
length=11
particles=3

total_crude=0
total_course=0

number_of_seeds=1
for seed in $(seq 1 $number_of_seeds)
do
  ./performance_test_crude_vs_coursegrain $sigma $length $seed $resolution $particles > out.txt
  value_crude=$(cat out.txt | grep "Crude Monte Carlo Run Time" | awk '{print $6}') 
  echo $value_crude
  total_crude=$(( $value_crude+$total_crude ))
  value_course=$(cat out.txt | grep "Course Monte Carlo Run Time" | awk '{print $6}') 
  echo $value_course
  total_course=$(( $value_course+$total_course )) 
done

echo $total_crude
echo $total_course

average_time_crude=$(bc <<< "scale = 6; ( $total_crude / $number_of_seeds )")
average_time_course=$(bc <<< "scale = 6; ( $total_course / $number_of_seeds )")

echo $average_time_crude
echo $average_time_course

echo $average_time_crude >> average_crude_time.txt
echo $average_time_course >> average_course_time.txt
