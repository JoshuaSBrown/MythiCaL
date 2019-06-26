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
		KMC_Site & createKMC_Site(const int & siteId);

		void reserve(size_t num_buckets) { sites_.reserve(num_buckets);}
    void addKMC_Site(KMC_Site& site);
    void addKMC_Sites(std::vector<KMC_Site>& sites);

    KMC_Site& getKMC_Site(const int & siteId);

//    std::unordered_map<int,KMC_Site> getKMC_Sites(std::vector<int> siteIds) const;
//    std::unordered_map<int,KMC_Site> getKMC_Sites() const;
    size_t size() const {return sites_.size();}

    void setClusterId(const int & siteId,const int &clusterId);

    int getClusterIdOfSite(const int & siteId) const;

    bool partOfCluster(const int & siteId) const;

    int getSmallestClusterId(std::vector<int> siteIds) const;

    bool exist(const int & siteId) const;

    bool isOccupied(const int & siteId) const;

    void vacate(const int & siteId);

    void occupy(const int & siteId);

    std::vector<int> getSiteIds() const; 

    double getDwellTime(const int &siteId);

    double getTimeConstant(const int &siteId) const;

    std::unordered_map<int,double> & getNeighborsAndRates(const int & siteId);

    Rate_Map getRates() const;

    double getFastestRateOffSite(const int & siteId) const;

    double getRateToNeighborOfSite(const int & siteId, const int & neighId) const;

    std::vector<int> getSiteIdsOfNeighbors(const int & siteId) const;
  private:
    std::unordered_map<int,KMC_Site> sites_;

};

}

#endif // KMCCOARSEGRAIN_KMC_SITE_CONTAINER_HPP
