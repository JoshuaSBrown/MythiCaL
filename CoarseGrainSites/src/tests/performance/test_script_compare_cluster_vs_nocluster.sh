#!/bin/bash

sigma=0.09
length=10
seed=11
walkers=3
threshold=1000
time=0.005
sample_rate=100


./performance_test_coarsegrain_cluster_vs_nocluster $sigma $length $seed $walkers $threshold $time $sample_rate> out.txt
value_coarse_nocluster=$(cat out.txt | grep "Coarse nocluster Monte Carlo Run Time" | awk '{print $7}') 
value_coarse=$(cat out.txt | grep "Coarse Monte Carlo Run Time" | awk '{print $6}') 

echo "No cluster time "$value_coarse_nocluster
echo "With cluster time "$value_coarse

if [ "$value_coarse" -lt "$value_coarse_nocluster" ]
then
  echo "Success"
  exit 0
else
  echo "Fail"
  exit 1
fi  
