
#include <unordered_map>

#include "kmc_cluster_container.hpp"

using namespace std;

namespace kmccoarsegrain {

  void KMC_Cluster_Container::addKMC_Cluster(KMC_Cluster& cluster) {
    if(exist(cluster.getId())){
      throw invalid_argument("Cannot add cluster as it is already stored in the"
          " container.");
    }
    clusters_.insert( std::unordered_map<int,KMC_Cluster>::value_type(cluster.getId(), KMC_Cluster(cluster))); 
  }

/*  void KMC_Cluster_Container::addKMC_Clusters(vector<KMC_Cluster>& clusters){
    for( auto & cluster : clusters){
      addKMC_Cluster(cluster);
    }
  }*/

  KMC_Cluster& KMC_Cluster_Container::getKMC_Cluster(const int & clusterId){
    if(exist(clusterId)==false){
      cerr << "Trying to access cluster with id " << clusterId << endl;
      throw invalid_argument("Cannot get cluster as it is not stored in the"
          " container.");
    }
    return clusters_.at(clusterId);
  }

  bool KMC_Cluster_Container::exist(const int & clusterId) const{
    if(clusters_.count(clusterId)){
      return true;
    }
    return false;
  }

  void KMC_Cluster_Container::erase(const int & clusterId){
    auto iter = clusters_.find(clusterId);
    if(iter!=clusters_.end()){
      clusters_.erase(iter);
    }
  }

  bool KMC_Cluster_Container::isOccupied(const int & clusterId) const{
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot determine if cluster is occupied as it"
          " is not stored in the container.");
    }
    return clusters_.at(clusterId).isOccupied();
  }

/*  void KMC_Cluster_Container::vacate(const int & clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot vacate cluster as it is not stored in the "
          "container.");
    }
    clusters_[clusterId].vacate();
  }*/

/*  void KMC_Cluster_Container::occupy(const int & clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot occupy cluster as it is not stored in the "
          "container.");
    }
    clusters_[clusterId].occupy();
  }*/

  vector<int> KMC_Cluster_Container::getClusterIds() const {
    vector<int> clusterids;
    for( const pair<int,KMC_Cluster> & cluster_pair : clusters_){
      clusterids.push_back(cluster_pair.first);
    }
    return clusterids;
  }

  double KMC_Cluster_Container::getDwellTime(const int & walker_id,const int & clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get cluster dwell time as it is not stored"
          " in the container.");
    }
    return clusters_.at(clusterId).getDwellTime(walker_id);
  }

/*  double KMC_Cluster_Container::getTimeConstant(const int & clusterId) const{
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get cluster time constant as it is not "
          "stored in the container.");
    }
    return clusters_.at(clusterId).getTimeConstant();
  }*/

  double KMC_Cluster_Container::getFastestRateOffCluster(const int & clusterId) const {
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get fastest rate off cluaster as it is not"
          " stored in the container.");
    }
    return clusters_.at(clusterId).getFastestRateOffCluster();
  }

  vector<int> KMC_Cluster_Container::getSiteIdsOfNeighbors(const int & clusterId) const {
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get site ids neighboring cluster as the "
          " cluster is not stored in the container.");
    }
    return clusters_.at(clusterId).getSiteIdsNeighboringCluster();
  }

  unordered_map<int,vector<int>> KMC_Cluster_Container::getSiteIdsOfClusters() const {
    unordered_map<int,vector<int>> clusters;
    for(const pair<int,KMC_Cluster> & cluster : clusters_){
      clusters[cluster.first]=cluster.second.getSiteIdsInCluster();
    }
    return clusters;
  }

  unordered_map<int,double> KMC_Cluster_Container::getResolutionOfClusters() const{
    unordered_map<int,double> clusters;
    for(const pair<int,KMC_Cluster> & cluster : clusters_){
      clusters[cluster.first]=cluster.second.getResolution();
    }
    return clusters;
  }

  unordered_map<int,double> KMC_Cluster_Container::getTimeIncrementOfClusters() const {
    unordered_map<int,double> clusters;
    for(const pair<int,KMC_Cluster> & cluster : clusters_){
      clusters[cluster.first]=cluster.second.getTimeIncrement();
    }
    return clusters;
  }

}

