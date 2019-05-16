
#include "kmc_site_container.hpp"

using namespace std;

namespace kmccoarsegrain {

  void KMC_Site_Container::addKMC_Site(KMC_Site& site){
    if(sites_.count(site.getId())){
      throw invalid_argument("Cannot add site it has already been added.");
    }
    sites_[site.getId()] = site;
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

  unordered_map<int,KMC_Site> KMC_Site_Container::getKMC_Sites(vector<int> siteIds){
    unordered_map<int,KMC_Site> sites;
    for( auto siteId : siteIds ){
      if(sites_.count(siteId)){
        sites[siteId] = sites_[siteId];
      }else{
        throw invalid_argument("Site is not found in the container.");
      }
    }
    return sites;
  }

  unordered_map<int,KMC_Site> KMC_Site_Container::getKMC_Sites(){
    return sites_;
  } 

  void KMC_Site_Container::setClusterId(int siteId, int clusterId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    sites_[siteId].setClusterId(clusterId);
  }

  int KMC_Site_Container::getClusterIdOfSite(int siteId) {
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    return sites_[siteId].getClusterId();
  }

  bool KMC_Site_Container::partOfCluster(int siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Site is not stored in the container.");
    }
    return sites_[siteId].partOfCluster();
  }

  int KMC_Site_Container::getSmallestClusterId(vector<int> siteIds){
    LOG("Getting the favored cluster Id", 1);
    int favoredClusterId = constants::unassignedId;
    for (auto siteId : siteIds) {
      int clusterId = sites_[siteId].getClusterId();
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
  
  bool KMC_Site_Container::isOccupied(const int & siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot determine if site is occupied as it is not"
          " stored in the container");
    }
    return sites_[siteId].isOccupied();
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
/*
  Rate_Map KMC_Site_Container::getInternalRates(vector<int> siteIds){
    
  }

  Rate_Map KMC_Site_Container::getExternalRates(vector<int> siteIds){

  }
*/
  vector<int> KMC_Site_Container::getSiteIds(){
    vector<int> siteIds;
    for(auto site : sites_ ){
      siteIds.push_back(site.first);
    }
    return siteIds;
  }
/*
  double KMC_Site_Container::getMaxTraverseTimeOfConnectedSites(vector<int> siteIds){
    auto edges = createEdges_(sites_, siteIds);
    auto nodes = createNodes_(siteIds);
    list<weak_ptr<Edge>> edges_weak(edges.begin(), edges.end());
    unordered_map<int, weak_ptr<GraphNode<string>>> nodes_weak;
    for (auto map_iter : nodes) nodes_weak[map_iter.first] = map_iter.second;

    auto graph_ptr =
      shared_ptr<Graph<string>>(new Graph<string>(edges_weak, nodes_weak));

    unordered_map<pair<int, int>, double,hash_functions::hash> verticesAndtimes =
      maxMinimumDistanceBetweenEveryVertex<string>(*graph_ptr);

    double maxtime = 0.0;
    for (auto verticesAndTime : verticesAndtimes) {
      if (verticesAndTime.second > maxtime) maxtime = verticesAndTime.second;
    }
    return maxtime;
  }
*/
  double KMC_Site_Container::getDwellTime(int siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get site dwell time as site is not in the "
          "container.");
    }
    return sites_[siteId].getDwellTime(constants::unassignedId);
  }

  double KMC_Site_Container::getTimeConstant(int siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get site time constant as site is not in "
          "the container.");
    }
    return sites_[siteId].getTimeConstant();
  }

  Rate_Map KMC_Site_Container::getRates(){
    Rate_Map rate_map;
    for( auto & site : sites_ ){
      rate_map[site.first] = site.second.getNeighborsAndRates();
    }
    return rate_map;
  }

  double KMC_Site_Container::getFastestRateOffSite(int siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get fastest rate off site as it is not "
          "stored in the container.");
    }
    return sites_[siteId].getFastestRate();
  }

  double KMC_Site_Container::getRateToNeighborOfSite(int siteId, int neighId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get rate from site to neighbor as site is "
          "not stored in the container");
    }
    return sites_[siteId].getRateToNeighbor(neighId);    
  }

  vector<int> KMC_Site_Container::getSiteIdsOfNeighbors(int siteId){
    if(sites_.count(siteId)==0){
      throw invalid_argument("Cannot get neighbor site ids from site as site is "
          "not stored in the container");
    }
    return sites_[siteId].getNeighborSiteIds();
  }
}
