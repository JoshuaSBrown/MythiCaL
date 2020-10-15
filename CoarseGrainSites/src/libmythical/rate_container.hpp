#ifndef MYTHICAL_RATE_CONTAINER_HPP
#define MYTHICAL_RATE_CONTAINER_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace mythical {

typedef std::unordered_map<int,std::unordered_map<int, double *>> Rate_Map;
/**
 * \brief Class is designed to store rates 
 **/
class Rate_Container {
  public:
    Rate_Container() {};
    Rate_Container(Rate_Map rates) : rates_(rates) {};

    void addRate(int siteId, int neighId,double * rate);
    void addRates(Rate_Map rates);
    double * getRate(int siteId, int neighId);

    size_t incomingRateCount(int siteId);
    size_t outgoingRateCount(int siteId);

    Rate_Map getIncomingRates(int siteId);
    Rate_Map getOutgoingRates(int siteId);

    std::vector<int> getSourceSiteIds(); 
    std::vector<int> getSinkSiteIds(); 

  private:
    Rate_Map rates_;
    std::unordered_set<int> getSiteIdsWithIncomingRates_();
    std::unordered_set<int> getSiteIdsWithOutgoingRates_();
};

}

#endif // MYTHICAL_RATE_CONTAINER_HPP
