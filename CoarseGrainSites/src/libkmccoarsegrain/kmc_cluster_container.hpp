#ifndef KMCCOARSEGRAIN_KMC_CLUSTER_CONTAINER_HPP
#define KMCCOARSEGRAIN_KMC_CLUSTER_CONTAINER_HPP

#include <unordered_map>

#include "log.hpp"
#include "kmc_rate_container.hpp"
#include "topologyfeatures/kmc_cluster.hpp"

namespace kmccoarsegrain {

/**
 * \brief Class is designed to store kmc sites
 *
 **/
class KMC_Cluster_Container {
  public:
    KMC_Cluster_Container() {};

    void addKMC_Cluster(KMC_Cluster& cluster);
    void addKMC_Clusters(std::vector<KMC_Cluster>& clusters);
    KMC_Cluster& getKMC_Cluster(int clusterId);

    size_t size();

    bool exist(int clusterId);
    void erase(int clusterId);
    bool isOccupied(int clusterId);
    void vacate(int clusterId);
    void occupy(int clusterId);

    std::vector<int> getClusterIds(); 
    double getDwellTime(int walker_id, int clusterId);
    double getTimeConstant(int clusterId);

    double getFastestRateOffCluster(int clusterId);
    std::vector<int> getSiteIdsOfNeighbors(int clusterId);

    std::vector<std::vector<int>> getSiteIdsOfClusters();

  private:
    std::unordered_map<int,KMC_Cluster> clusters_;

};

}

#endif // KMCCOARSEGRAIN_CLUSTER_SITE_CONTAINER_HPP
