#ifndef KMCCOARSEGRAIN_KMC_CLUSTER_CONTAINER_HPP
#define KMCCOARSEGRAIN_KMC_CLUSTER_CONTAINER_HPP

#include <unordered_map>
#include <vector>

#include "log.hpp"
#include "kmc_cluster.hpp"

namespace kmccoarsegrain {

/**
 * \brief Class is designed to store kmc sites
 *
 **/
class KMC_Cluster_Container {
  public:
    KMC_Cluster_Container() {};

    void addKMC_Cluster(KMC_Cluster& cluster);

    KMC_Cluster& getKMC_Cluster(const int & clusterId);

    size_t size() const { return clusters_.size();}

    bool exist(const int & clusterId) const;

    void erase(const int & clusterId);

    bool isOccupied(const int & clusterId) const;

   // void vacate(const int & clusterId);

   // void occupy(const int & clusterId);

    std::vector<int> getClusterIds() const; 
    double getDwellTime(const int & walker_id, const int & clusterId);

    //double getTimeConstant(const int & clusterId) const;

    double getFastestRateOffCluster(const int & clusterId) const;
    std::vector<int> getSiteIdsOfNeighbors(const int & clusterId) const;

    std::unordered_map<int,double> getResolutionOfClusters() const;
    std::unordered_map<int,double> getTimeIncrementOfClusters() const;
    std::unordered_map<int,std::vector<int>> getSiteIdsOfClusters() const;

  private:
    std::unordered_map<int,KMC_Cluster> clusters_;

};

}

#endif // KMCCOARSEGRAIN_CLUSTER_SITE_CONTAINER_HPP
