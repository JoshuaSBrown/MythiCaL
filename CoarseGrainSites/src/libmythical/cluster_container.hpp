#ifndef MYTHICAL_CLUSTER_CONTAINER_HPP
#define MYTHICAL_CLUSTER_CONTAINER_HPP

#include <unordered_map>
#include <vector>

#include "log.hpp"
#include "rate_container.hpp"
#include "topologyfeatures/cluster.hpp"

namespace mythical {

/**
 * \brief Class is designed to store kmc sites
 *
 **/
class Cluster_Container {
  public:
    Cluster_Container() {};

    void addCluster(Cluster& cluster);
    void addClusters(std::vector<Cluster>& clusters);
    Cluster& getCluster(int clusterId);

    size_t size() const { return clusters_.size();}

    bool exist(const int & clusterId) const;
    void erase(int clusterId);
    bool isOccupied(const int & clusterId);
    void vacate(const int & clusterId);
    void occupy(const int & clusterId);

    std::vector<int> getClusterIds(); 
    double getDwellTime(int walker_id, int clusterId);
    double getTimeConstant(int clusterId);

    double getFastestRateOffCluster(int clusterId);
    std::vector<int> getSiteIdsOfNeighbors(int clusterId);

    std::unordered_map<int,double> getResolutionOfClusters();
    std::unordered_map<int,double> getTimeIncrementOfClusters();
    std::unordered_map<int,std::vector<int>> getSiteIdsOfClusters();

  private:
    std::unordered_map<int,Cluster> clusters_;

};

}

#endif // MYTHICAL_CLUSTER_CONTAINER_HPP
