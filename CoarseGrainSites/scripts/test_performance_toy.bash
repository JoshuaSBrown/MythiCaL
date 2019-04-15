#!/bin/bash

rm out.txt

charges=100000
threshold=10
time1=1
echo "time 2 | For time increment = 1                                       | For time increment = 5                                       | For time increment = 50" > performance2.txt
echo "       | CMC time | CME time | CME cluster Res | CME cluster Time inc | CMC time | CME time | CME cluster Res | CME cluster Time inc | CMC time | CME time | CME cluster Res | CME cluster Time inc |" >> performance2.txt
for time2 in {1..100..2}
do
  ./performance_test_coarsegrain_cluster_vs_nocluster_toy $time1 $time2 $charges $threshold 1 > out.txt
  crude_time=$(grep "Crude Monte Carlo Run Time:" out.txt | awk '{ print $6}')
  coarse_time_res_inc=$(grep "Coarse Monte Carlo Run Time:" out.txt | awk '{print $6"   "$9"   "$12}') 
  ./performance_test_coarsegrain_cluster_vs_nocluster_toy $time1 $time2 $charges $threshold 5 > out2.txt
  crude_time2=$(grep "Crude Monte Carlo Run Time:" out2.txt | awk '{ print $6}')
  coarse_time_res_inc2=$(grep "Coarse Monte Carlo Run Time:" out2.txt | awk '{print $6"   "$9"   "$12}') 
  ./performance_test_coarsegrain_cluster_vs_nocluster_toy $time1 $time2 $charges $threshold 50 > out3.txt
  crude_time3=$(grep "Crude Monte Carlo Run Time:" out3.txt | awk '{ print $6}')
  coarse_time_res_inc3=$(grep "Coarse Monte Carlo Run Time:" out3.txt | awk '{print $6"   "$9"   "$12}') 

  echo "$time2   $crude_time   $coarse_time_res_inc   $crude_time2   $coarse_time_res_inc2   $crude_time3   $coarse_time_res_inc3" >> performance2.txt 
done
