#ifndef KMCCOURSEGRAIN_KMC_CLUSTER_CONTAINER_HPP
#define KMCCOURSEGRAIN_KMC_CLUSTER_CONTAINER_HPP

#include <unordered_map>

#include "log.hpp"
#include "kmc_rate_container.hpp"
#include "topologyfeatures/kmc_cluster.hpp"

namespace kmccoursegrain {

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

    bool isOccupied(int clusterId);
    void vacate(int clusterId);
    void occupy(int clusterId);

    std::vector<int> getClusterIds(); 
    double getDwellTime(int clusterId);
    double getTimeConstant(int clusterId);

    double getFastestRateOffCluster(int clusterId);
    std::vector<int> getSiteIdsOfNeighbors(int clusterId);
  private:
    std::unordered_map<int,KMC_Cluster> clusters_;

};

}

#endif // KMCCOURSEGRAIN_CLUSTER_SITE_CONTAINER_HPP
