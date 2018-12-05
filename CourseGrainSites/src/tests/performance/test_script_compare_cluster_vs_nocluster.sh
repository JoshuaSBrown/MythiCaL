#!/bin/bash

sigma=0.09
length=10
seed=11
walkers=3
threshold=1000
time=0.005
sample_rate=100


./performance_test_coursegrain_cluster_vs_nocluster $sigma $length $seed $walkers $threshold $time $sample_rate> out.txt
value_course_nocluster=$(cat out.txt | grep "Course nocluster Monte Carlo Run Time" | awk '{print $7}') 
value_course=$(cat out.txt | grep "Course Monte Carlo Run Time" | awk '{print $6}') 

echo "No cluster time "$value_course_nocluster
echo "With cluster time "$value_course

if [ "$value_course" -lt "$value_course_nocluster" ]
then
  echo "Success"
  exit 0
else
  echo "Fail"
  exit 1
fi  
