#ifndef KMCCOURSEGRAIN_KMC_SITE_CONTAINER_HPP
#define KMCCOURSEGRAIN_KMC_SITE_CONTAINER_HPP

#include <unordered_map>

#include "log.hpp"
#include "kmc_rate_container.hpp"
#include "topologyfeatures/kmc_site.hpp"

namespace kmccoursegrain {

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

    /**
     * \brief get the rates to the sites stored in the container 
     **/
//    Rate_Map getInternalRates(std::vector<int> siteIds);
    /**
     * \brief get the rates to the sites that are not stored in the container
     **/
//    Rate_Map getExternalRates(std::vector<int> siteIds);

    std::vector<int> getSiteIds(); 
 //   double getMaxTraverseTimeOfConnectedSites(std::vector<int> siteIds);
    double getDwellTime(int siteId);
    double getTimeConstant(int siteId);

    Rate_Map getRates();

    double getFastestRateOffSite(int siteId);
    double getRateToNeighborOfSite(int siteId, int neighId);
    std::vector<int> getSiteIdsOfNeighbors(int siteId);
//    double getProbabilityOfHoppingToNeighbor(int siteId, int neighId);
//    void moveSitesFrom(KMC_Site_Container& container);
  private:
    std::unordered_map<int,KMC_Site> sites_;

};

}

#endif // KMCCOURSEGRAIN_KMC_SITE_CONTAINER_HPP
