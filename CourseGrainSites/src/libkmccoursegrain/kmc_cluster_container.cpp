
#include <unordered_map>

#include "kmc_cluster_container.hpp"

using namespace std;

namespace kmccoursegrain {

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
      throw invalid_argument("Cannot get cluster as it is not stored in the"
          " container.");
    }
    return clusters_[clusterId];
  }

  size_t KMC_Cluster_Container::size() { return clusters_.size(); }

  bool KMC_Cluster_Container::exist(int clusterId){
    if(clusters_.count(clusterId)){
      return true;
    }
    return false;
  }

  bool KMC_Cluster_Container::isOccupied(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot determine if cluster is occupied as it"
          " is not stored in the container.");
    }
    return clusters_[clusterId].isOccupied();
  }

  void KMC_Cluster_Container::vacate(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot vacate cluster as it is not stored in the "
          "container.");
    }
    clusters_[clusterId].vacate();
  }

  void KMC_Cluster_Container::occupy(int clusterId){
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

  double KMC_Cluster_Container::getDwellTime(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get cluster dwell time as it is not stored"
          " in the container.");
    }
    return clusters_[clusterId].getDwellTime();
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

  std::vector<int> KMC_Cluster_Container::getSiteIdsOfNeighbors(int clusterId){
    if(exist(clusterId)==false){
      throw invalid_argument("Cannot get site ids neighboring cluster as the "
          " cluster is not stored in the container.");
    }
    return clusters_[clusterId].getSiteIdsNeighboringCluster();
  }


}
