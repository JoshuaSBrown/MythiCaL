#ifndef KMCCOARSEGRAIN_KMC_DYNAMIC_TOPOLOGY_HPP
#define KMCCOARSEGRAIN_KMC_DYNAMIC_TOPOLOGY_HPP

#include <unordered_map>
#include <vector>

#include "log.hpp"
#include "kmc_rate_container.hpp"
#include "kmc_site_container.hpp"
#include "kmc_cluster_container.hpp"
#include "topologyfeatures/kmc_site.hpp"

namespace kmccoarsegrain {

class KMC_Dynamic_Topology;
class KMC_Cluster;

KMC_TopologyFeature * returnSite(KMC_Dynamic_Topology & dyn_top,int siteId);
KMC_TopologyFeature * returnCluster(KMC_Dynamic_Topology & dyn_top,int siteId);
KMC_TopologyFeature * createSite(KMC_Dynamic_Topology & dyn_top,int siteId);
/**
 * \brief Class is in charge of the topology objects of the simulation
 *
 * The topology objects include clusters and sites. The dynamic topology is 
 * designed such that memory is only allocated to a topology object when it is
 * used. In essence, it assumes not all sites will be utilized so memory is
 * not allocated to them until they are used. 
 *
 **/
class KMC_Dynamic_Topology {
  public:
    KMC_Dynamic_Topology() : seed_set_(false), seed_(0) {};

    void setRates(std::unordered_map<int,std::unordered_map<int,double>> & rates);

    void setRandomSeed(const unsigned long seed);

    int getClusterIdOfSite(int siteId);

    int getVisitFrequencyOfSite(int siteId);

    void addKMC_Cluster(KMC_Cluster & cluster);

    KMC_Site & getKMC_Site(const int siteId);
    KMC_Cluster & getKMC_Cluster(const int clusterId);

    int getFavoredClusterId(std::vector<int> siteIds);

    bool siteExist(int siteId);

    void setSitesClusterId(int siteId, int clusterId);

    void mergeClusters(int clusterId1, int clusterId2);

    double getFastestRateOffSite(const int siteId);
    double getFastestRateOffCluster(const int clusterId);

    std::unordered_map<int,double> getResolutionOfClusters();    

    std::unordered_map<int,double> & getSiteRates(const int siteId);    
    double getRate(const int siteId1, const int siteId2);

    bool partOfCluster(const int siteId);

    std::unordered_map<int,std::vector<int>> getClusters();

    std::unordered_map<int,int> getClustersOfSites(const std::vector<int> & siteIds);
    typedef KMC_TopologyFeature * (*top_feature)(KMC_Dynamic_Topology & top, int siteId);

    struct DefaultTopologyFeatureFunction {
      top_feature feature = &createSite;
    };

    std::unordered_map<int,double> getTimeIncrementOfClusters();

    double getTimeConstantFromSitesToNeighbors(const std::vector<int> & siteIds) const;
    std::unordered_map<int,DefaultTopologyFeatureFunction> features;
  private:

    bool seed_set_;
    unsigned long seed_;
    std::unordered_map<int,std::unordered_map<int,double>> * rates_;
    //friend class BasinExplorer;

    friend KMC_TopologyFeature * returnSite(KMC_Dynamic_Topology & dyn_top,int siteId);
    friend KMC_TopologyFeature * returnCluster(KMC_Dynamic_Topology & dyn_top,int siteId);
    friend KMC_TopologyFeature * createSite(KMC_Dynamic_Topology & dyn_top,int siteId);

    /// Stores smart pointers to all the sites
    KMC_Site_Container sites_;

    /// Stores smart pointers to all the clusters
    KMC_Cluster_Container clusters_;
};

}

#endif // KMCCOARSEGRAIN_KMC_DYNAMIC_TOPOLOGY_HPP
