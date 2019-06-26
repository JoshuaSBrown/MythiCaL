#ifndef KMCCOARSEGRAIN_KMC_RATE_CONTAINER_HPP
#define KMCCOARSEGRAIN_KMC_RATE_CONTAINER_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace kmccoarsegrain {

typedef std::unordered_map<int,std::unordered_map<int,double> *> Rate_Map;
/**
 * \brief Class is designed to store rates 
 **/
class KMC_Rate_Container {
  public:
    KMC_Rate_Container() {};
    KMC_Rate_Container(Rate_Map rates) : rates_(rates) {};

    //void addRate(int siteId, int neighId,double * rate);
    void addRates(Rate_Map rates);
    double getRate(int siteId, int neighId);

    size_t incomingRateCount(int siteId) const;
    size_t outgoingRateCount(int siteId) const;

    //Rate_Map getIncomingRates(int siteId);
    //Rate_Map getOutgoingRates(int siteId);

    std::vector<int> getSourceSiteIds(); 
    std::vector<int> getSinkSiteIds(); 

  private:
    Rate_Map rates_;
    std::unordered_set<int> getSiteIdsWithIncomingRates_();
    std::unordered_set<int> getSiteIdsWithOutgoingRates_();
};

}

#endif // KMCCOARSEGRAIN_KMC_RATE_CONTAINER_HPP
