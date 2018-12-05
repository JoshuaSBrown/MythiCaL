sigma=0.07
size=30
seed=1
threshold=500
particles=30
time=0.0001
field=1.0
sample_rate=100

./regression_test_ToF $sigma $size $seed $particles $threshold $time $sample_rate $field

# Test when clusters are used
./regression_test_1Dsystem 8000 1 10 0.002 1000 10 > out.txt

# Test when clusters are essentially turned off
./regression_test_1Dsystem 8000 1 10000000 0.002 1000 10 > out2.txt

cluster_deviation=$(cat out.txt | grep "Standard deviation" | awk '{print $3}')
cluster_mean=$(cat out.txt | grep "Mean" | awk '{print $2}')
number_clusters=$(cat out.txt | grep "Number of clusters" | awk '{print $5}')

nocluster_deviation=$(cat out2.txt | grep "Standard deviation" | awk '{print $3}')
nocluster_mean=$(cat out2.txt | grep "Mean" | awk '{print $2}')

# Ensure that clusters were found
if [ "$cluster_mean" -lt "0" ] 
then
  echo "No clusters were found when the algorithm should have detected some"
  exit 1
fi

# Ensure that with cluters the signal is within 10% of the algorithm when no
# clusters are detected

upper_range_dev=$(echo "$nocluster_deviation*1.1" | bc)
lower_range_dev=$(echo "$nocluster_deviation*0.9" | bc)

upper_range_mean=$(echo "$nocluster_mean*1.1" | bc)
lower_range_mean=$(echo "$nocluster_mean*0.9" | bc)

if [ "$cluster_deviation" -gt "$upper_range_dev" ] 
then
  echo "Cluster noise deviation is too large when compared to algorithm without"
  echo "clusters."
  exit 1
fi    

if [ "$cluster_mean" -gt "$upper_range_mean" ] 
then
  echo "Cluster noise mean is too large when compared to algorithm without"
  echo "clusters."
  exit 1
fi    

if [ "$cluster_deviation" -lt "$lower_range_dev" ]
then
  echo "Cluster noise deviation is too low when compared to algorithm without "
  echo "clusters."
  exit 1
fi  

if [ "$cluster_mean" -lt "$lower_range_mean" ]
then
  echo "Cluster noise mean is too low when compared to algorithm without "
  echo "clusters."
  exit 1
fi  

exit 0
