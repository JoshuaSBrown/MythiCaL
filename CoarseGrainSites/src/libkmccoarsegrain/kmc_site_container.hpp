#ifndef KMCCOARSEGRAIN_KMC_SITE_CONTAINER_HPP
#define KMCCOARSEGRAIN_KMC_SITE_CONTAINER_HPP

#include <unordered_map>

#include "log.hpp"
#include "kmc_rate_container.hpp"
#include "topologyfeatures/kmc_site.hpp"

namespace kmccoarsegrain {

/**
 * \brief Class is designed to store kmc sites
 *
 **/
class KMC_Site_Container {
  public:
    KMC_Site_Container() {};

    void addKMC_Site(KMC_Site& site);
    void addKMC_Sites(std::vector<KMC_Site>& sites);
    KMC_Site& getKMC_Site(int siteId);

    std::unordered_map<int,KMC_Site> getKMC_Sites(std::vector<int> siteIds);
    std::unordered_map<int,KMC_Site> getKMC_Sites();
    size_t size();

    void setClusterId(int siteId, int clusterId);
    int getClusterIdOfSite(int siteId);
    bool partOfCluster(int siteId);
    int getSmallestClusterId(std::vector<int> siteIds);

    bool exist(int siteId);

    bool isOccupied(int siteId);
    void vacate(int siteId);
    void occupy(int siteId);

    std::vector<int> getSiteIds(); 
    double getDwellTime(int siteId);
    double getTimeConstant(int siteId);

    Rate_Map getRates();

    double getFastestRateOffSite(int siteId);
    double getRateToNeighborOfSite(int siteId, int neighId);
    std::vector<int> getSiteIdsOfNeighbors(int siteId);
  private:
    std::unordered_map<int,KMC_Site> sites_;

};

}

#endif // KMCCOARSEGRAIN_KMC_SITE_CONTAINER_HPP
