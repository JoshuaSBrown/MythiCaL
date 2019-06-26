#include <utility>
#include "kmc_site_container.hpp"

using namespace std;

namespace kmccoarsegrain {

  void KMC_Site_Container::addKMC_Site(KMC_Site& site){
    if(sites_.count(site.getId())){
      throw invalid_argument("Cannot add site it has already been added.");
    }
    sites_[site.getId()] = site;
  }

	KMC_Site &  KMC_Site_Container::createKMC_Site(const int & siteId){
    if(sites_.count(siteId)){
      throw invalid_argument("Cannot add site it has already been added.");
    }
		KMC_Site & site = sites_[siteId];
		site.setId(siteId);
		return site;
  }

  void KMC_Site_Container::addKMC_Sites(vector<KMC_Site>& sites){
    for ( KMC_Site & site : sites ){
      addKMC_Site(site);
    }
  }

  KMC_Site& KMC_Site_Container::getKMC_Site(const int & siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    return sites_[siteId];
  }

  void KMC_Site_Container::setClusterId(const int & siteId, const int & clusterId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    sites_[siteId].setClusterId(clusterId);
  }

  int KMC_Site_Container::getClusterIdOfSite(const int & siteId) const{
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    return sites_.at(siteId).getClusterId();
  }

  bool KMC_Site_Container::partOfCluster(const int & siteId) const{
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    return sites_.at(siteId).partOfCluster();
  }

  int KMC_Site_Container::getSmallestClusterId(vector<int> siteIds) const {
    LOG("Getting the favored cluster Id", 1);
    int favoredClusterId = constants::unassignedId;
    for (const int & siteId : siteIds) {
      int clusterId = sites_.at(siteId).getClusterId();
      if (favoredClusterId == constants::unassignedId) {
        favoredClusterId = clusterId;
      } else if (clusterId != constants::unassignedId &&
          clusterId < favoredClusterId) {
        favoredClusterId = clusterId;
      }
    }
    return favoredClusterId;
  }

  bool KMC_Site_Container::exist(const int & siteId) const{
    return sites_.count(siteId)!=0;
  }
  
  bool KMC_Site_Container::isOccupied(const int & siteId) const{
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot determine if site is occupied as it is not"
          " stored in the container");
    }
    return sites_.at(siteId).isOccupied();
  }

  void KMC_Site_Container::vacate(const int & siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot vacate as site is not stored in the "
          "container.");
    }
    sites_[siteId].vacate();
  }

  void KMC_Site_Container::occupy(const int & siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot occupy site as it is not stored in the "
          "container.");
    }
    sites_[siteId].occupy();
  }

  vector<int> KMC_Site_Container::getSiteIds() const{
    vector<int> siteIds;
    for(const pair<int,KMC_Site> & site : sites_ ){
      siteIds.push_back(site.first);
    }
    return siteIds;
  }

  double KMC_Site_Container::getDwellTime(const int & siteId) {
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get site dwell time as site is not in the "
          "container.");
    }
    return sites_.at(siteId).getDwellTime(constants::unassignedId);
  }

  double KMC_Site_Container::getTimeConstant(const int & siteId) const {
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get site time constant as site is not in "
          "the container.");
    }
    return sites_.at(siteId).getTimeConstant();
  }

  Rate_Map KMC_Site_Container::getRates() const {
    Rate_Map rate_map;
    for( const pair<int,KMC_Site> & site : sites_ ){
      rate_map[site.first] = &(site.second.getNeighborsAndRates());
    }
    return rate_map;
  }

  std::unordered_map<int,double> & KMC_Site_Container::getNeighborsAndRates(const int & siteId){
    return sites_[siteId].getNeighborsAndRates();
  }

  double KMC_Site_Container::getFastestRateOffSite(const int & siteId) const{
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get fastest rate off site as it is not "
          "stored in the container.");
    }
    return sites_.at(siteId).getFastestRate();
  }

  double KMC_Site_Container::getRateToNeighborOfSite(const int & siteId, const int & neighId) const{
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get rate from site to neighbor as site is "
          "not stored in the container");
    }
    return sites_.at(siteId).getRateToNeighbor(neighId);    
  }

  vector<int> KMC_Site_Container::getSiteIdsOfNeighbors(const int & siteId) const{
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get neighbor site ids from site as site is "
          "not stored in the container");
    }
    return sites_.at(siteId).getNeighborSiteIds();
  }
}
