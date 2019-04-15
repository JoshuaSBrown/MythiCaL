
#include <unordered_map>

#include "kmc_cluster_container.hpp"

using namespace std;

namespace kmccoarsegrain {

  void KMC_Cluster_Container::addKMC_Cluster(KMC_Cluster& cluster) {
    if(exist(cluster.getId())){
      throw invalid_argument("Cannot add cluster as it is already stored in the"
          " container.");
    }
    clusters_[cluster.getId()]=cluster; 
  }

  void KMC_Cluster_Container::addKMC_Clusters(vector<KMC_Cluster>& clusters){
    for( auto & cluster : clusters){
      addKMC_Cluster(cluster);
    }
  }

  KMC_Cluster& KMC_Cluster_Container::getKMC_Cluster(int clusterId){
    if(exist(clusterId)==false){
      cerr << "Trying to access cluster with id " << clusterId << endl;
      throw invalid_argument("Cannot get cluster as it is not stored in the"
          " container.");
    }
    return clusters_[clusterId];
  }

  bool KMC_Cluster_Container::exist(const int & clusterId) const{
    if(clusters_.count(clusterId)){
      return true;
    }
    return false;
  }

  void KMC_Cluster_Container::erase(int clusterId){
    auto iter = clusters_.find(clusterId);
    if(iter!=clusters_.end()){
      clusters_.erase(iter);
    }
  }

  bool KMC_Cluster_Container::isOccupied(const int & clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot determine if cluster is occupied as it"
          " is not stored in the container.");
    }
    return clusters_[clusterId].isOccupied();
  }

  void KMC_Cluster_Container::vacate(const int & clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot vacate cluster as it is not stored in the "
          "container.");
    }
    clusters_[clusterId].vacate();
  }

  void KMC_Cluster_Container::occupy(const int & clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot occupy cluster as it is not stored in the "
          "container.");
    }
    clusters_[clusterId].occupy();
  }

  vector<int> KMC_Cluster_Container::getClusterIds(){
    vector<int> clusterids;
    for( auto cluster_pair : clusters_){
      clusterids.push_back(cluster_pair.first);
    }
    return clusterids;
  }

  double KMC_Cluster_Container::getDwellTime(int walker_id, int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get cluster dwell time as it is not stored"
          " in the container.");
    }
    return clusters_[clusterId].getDwellTime(walker_id);
  }

  double KMC_Cluster_Container::getTimeConstant(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get cluster time constant as it is not "
          "stored in the container.");
    }
    return clusters_[clusterId].getTimeConstant();
  }

  double KMC_Cluster_Container::getFastestRateOffCluster(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get fastest rate off cluaster as it is not"
          " stored in the container.");
    }
    return clusters_[clusterId].getFastestRateOffCluster();
  }

  vector<int> KMC_Cluster_Container::getSiteIdsOfNeighbors(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get site ids neighboring cluster as the "
          " cluster is not stored in the container.");
    }
    return clusters_[clusterId].getSiteIdsNeighboringCluster();
  }

  unordered_map<int,vector<int>> KMC_Cluster_Container::getSiteIdsOfClusters(){
    unordered_map<int,vector<int>> clusters;
    for(auto cluster : clusters_){
      clusters[cluster.first]=cluster.second.getSiteIdsInCluster();
    }
    return clusters;
  }

  unordered_map<int,int> KMC_Cluster_Container::getResolutionOfClusters(){
    unordered_map<int,int> clusters;
    for(auto cluster : clusters_){
      clusters[cluster.first]=cluster.second.getResolution();
    }
    return clusters;
  }

  unordered_map<int,double> KMC_Cluster_Container::getTimeIncrementOfClusters(){
    unordered_map<int,double> clusters;
    for(auto cluster : clusters_){
      clusters[cluster.first]=cluster.second.getTimeIncrement();
    }
    return clusters;
  }

}

