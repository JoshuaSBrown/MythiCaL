#ifndef MYTHICAL_SITE_CONTAINER_HPP
#define MYTHICAL_SITE_CONTAINER_HPP

#include <unordered_map>

#include "log.hpp"
#include "rate_container.hpp"
#include "topologyfeatures/site.hpp"

namespace mythical {

/**
 * \brief Class is designed to store kmc sites
 *
 **/
class Site_Container {
  public:
    Site_Container() {};

    void addSite(Site& site);
    void addSites(std::vector<Site>& sites);
    Site& getSite(const int & siteId);

    std::unordered_map<int,Site> getSites(std::vector<int> siteIds);
    std::unordered_map<int,Site> getSites();
    size_t size() const {return sites_.size();}

    void setClusterId(int siteId, int clusterId);
    int getClusterIdOfSite(int siteId);
    bool partOfCluster(int siteId);
    int getSmallestClusterId(std::vector<int> siteIds);

    bool exist(const int & siteId) const;

    bool isOccupied(const int & siteId);
    void vacate(const int & siteId);
    void occupy(const int & siteId);

    std::vector<int> getSiteIds(); 
    double getDwellTime(int siteId);
    double getTimeConstant(int siteId);

    Rate_Map getRates();

    double getFastestRateOffSite(int siteId);
    double getRateToNeighborOfSite(int siteId, int neighId);
    std::vector<int> getSiteIdsOfNeighbors(int siteId);
  private:
    std::unordered_map<int,Site> sites_;

};

}

#endif // MYTHICAL_SITE_CONTAINER_HPP
