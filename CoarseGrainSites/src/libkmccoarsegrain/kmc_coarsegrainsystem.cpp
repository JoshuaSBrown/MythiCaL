
#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_set>

#include "../../include/kmccoarsegrain/kmc_coarsegrainsystem.hpp"
#include "../../include/kmccoarsegrain/kmc_constants.hpp"
#include "../../include/kmccoarsegrain/kmc_walker.hpp"

#include "topologyfeatures/kmc_topology_feature.hpp"
#include "topologyfeatures/kmc_cluster.hpp"
#include "topologyfeatures/kmc_site.hpp"
#include "log.hpp"
#include "kmc_basin_explorer.hpp"
#include "kmc_graph_library_adapter.hpp"
#include "kmc_site_container.hpp"
#include "kmc_cluster_container.hpp"

#include "../../../UGLY/include/ugly/pair_hash.hpp"
#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph.hpp"
#include "../../../UGLY/include/ugly/graph_algorithms.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"
#include "../../../UGLY/include/ugly/graphvisitor/graphvisitor_smallest_known_value.hpp"

using namespace std;
using namespace std::chrono;
using namespace ugly;
using namespace ugly::graphalgorithms;

namespace kmccoarsegrain {
  /****************************************************************************
   * Private Internal Function Declarations
   ****************************************************************************/

  bool compare(const pair<int,double> &x, const pair<int,double> &y){
    return x.second>y.second;
  }

  unordered_map<int, shared_ptr<GraphNode<string>>> createNode_(int siteIds);

  size_t countUniqueClusters(const unordered_map<int,int> & sites_and_clusters);
  int getFavoredClusterId(unordered_map<int,int> sites_and_clusters);

  /****************************************************************************
   * Public Facing Functions
   ****************************************************************************/

  KMC_CoarseGrainSystem::KMC_CoarseGrainSystem() :
    performance_ratio_(1.00),
    seed_set_(false),
    seed_(0),
    time_resolution_set_(false),
    minimum_coarse_graining_resolution_(2),
    iteration_(0),
    iteration_threshold_(1000),
    iteration_threshold_min_(1000){
      sites_ = unique_ptr<KMC_Site_Container>( new KMC_Site_Container );
      clusters_ = unique_ptr<KMC_Cluster_Container>( new KMC_Cluster_Container );
    }

  KMC_CoarseGrainSystem::~KMC_CoarseGrainSystem(){
  }

  double KMC_CoarseGrainSystem::getTimeResolution() { 
    if(!time_resolution_set_){
      throw runtime_error("Cannot get the time resolution as it has not yet "
          "been set.");
    }
    return time_resolution_; 
  }

  void KMC_CoarseGrainSystem::setTimeResolution(double time_resolution){
    if(time_resolution<=0.0){
      throw invalid_argument("The time resolution must be a positive value.");
    }
    time_resolution_set_ = true;
    time_resolution_ = time_resolution;
  }

  void KMC_CoarseGrainSystem::initializeSystem(unordered_map<int, unordered_map<int, double>>& ratesOfAllSites) {

    LOG("Initializeing system", 1);

    if(!time_resolution_set_){
      throw runtime_error("You must first set the time resolution of the system "
          "before you can initialize the system.");
    }

    for (auto it = ratesOfAllSites.begin(); it != ratesOfAllSites.end(); ++it) {
      KMC_Site site;
      site.setId(it->first);

      site.setRatesToNeighbors(it->second);
      if (seed_set_) {
        site.setRandomSeed(seed_);
        ++seed_;
      }
      sites_->addKMC_Site(site);
      topology_features_[it->first] = &(sites_->getKMC_Site(it->first));
    }

    // Address sites that will act as drains with no rates off of them
    unordered_set<int> drain_sites;
    for (pair<const int,unordered_map<int,double>> & sites_and_rates : ratesOfAllSites){
      for(const pair<int,double> & site_and_rate : sites_and_rates.second ){
        if(ratesOfAllSites.count(site_and_rate.first)==0){
          drain_sites.insert(site_and_rate.first);
        }
      }
    }

    for( const int & drain_site_id : drain_sites ){
      KMC_Site site;
      site.setId(drain_site_id);
      sites_->addKMC_Site(site);
      topology_features_[drain_site_id] = &(sites_->getKMC_Site(drain_site_id));
    }
  }

  int KMC_CoarseGrainSystem::getVisitFrequencyOfSite(int siteId){
    if(sites_->exist(siteId)==false){
      throw invalid_argument("Site is not stored in the coarse grained system you"
          " cannot retrieve it's visit frequency.");
    }

    int visits = sites_->getKMC_Site(siteId).getVisitFrequency();
    if(sites_->partOfCluster(siteId)){
      int cluster_id = sites_->getClusterIdOfSite(siteId);
      visits += clusters_->getKMC_Cluster(cluster_id).getVisitFrequency(siteId);
    }
    return visits;
  }

  void KMC_CoarseGrainSystem::initializeWalkers(vector<pair<int,KMC_Walker>>& walkers) {

    LOG("Initializeing walkers", 1);

    if (topology_features_.size() == 0) {
      throw runtime_error(
          "You must first initialize the system before you "
          "can initialize the walkers");
    }

    for ( size_t index = 0; index<walkers.size(); ++index){
      int siteId = walkers.at(index).second.getIdOfSiteCurrentlyOccupying();
      if (siteId == constants::unassignedId) {
        throw runtime_error(
            "You must first place the walker on a known site"
            " before the walker can be initialized.");
      }
      topology_features_[siteId]->occupy();

      auto hopTime = topology_features_[siteId]->getDwellTime(walkers.at(index).first);
      int newId = topology_features_[siteId]->pickNewSiteId(walkers.at(index).first);
      walkers.at(index).second.setDwellTime(hopTime);
      walkers.at(index).second.setPotentialSite(newId);
    }
  }

  void KMC_CoarseGrainSystem::setMinCoarseGrainIterationThreshold(int threshold_min) {
    LOG("Setting minimum coarse graining threshold", 1);
    iteration_threshold_min_ = threshold_min;
    iteration_threshold_ = threshold_min;
  }

  void KMC_CoarseGrainSystem::setRandomSeed(const unsigned long seed) {
    if (topology_features_.size() != 0) {
      throw runtime_error(
          "For the random seed to have an affect, it must be "
          "set before initializeSystem is called");
    }
    seed_ = seed;
    seed_set_ = true;
  }

  void KMC_CoarseGrainSystem::removeWalkerFromSystem(pair<int,KMC_Walker>& walker) {
    removeWalkerFromSystem(walker.first,walker.second);
  }

  void KMC_CoarseGrainSystem::removeWalkerFromSystem(int & walker_id, KMC_Walker& walker) {
    LOG("Walker is being removed from system", 1);
    auto siteId = walker.getIdOfSiteCurrentlyOccupying();
    topology_features_[siteId]->removeWalker(walker_id,siteId);
  }

  int KMC_CoarseGrainSystem::getClusterIdOfSite(int siteId) {
    return sites_->getClusterIdOfSite(siteId);
  }

  void KMC_CoarseGrainSystem::hop(pair<const int,KMC_Walker>& walker) {
    hop(walker.first,walker.second);
  }

  void KMC_CoarseGrainSystem::hop(const int & walker_id, KMC_Walker & walker) {
    const int & siteId = walker.getIdOfSiteCurrentlyOccupying();
    const int & siteToHopToId = walker.getPotentialSite();
    KMC_TopologyFeature *& feature = topology_features_[siteId];
    KMC_TopologyFeature *& feature_to_hop_to = topology_features_[siteToHopToId];

    if(!feature_to_hop_to->isOccupied(siteToHopToId)){
      feature->vacate(siteId);
      feature_to_hop_to->occupy(siteToHopToId);

      walker.occupySite(siteToHopToId);
      walker.setDwellTime(feature_to_hop_to->getDwellTime(walker_id));
      walker.setPotentialSite(feature_to_hop_to->pickNewSiteId(walker_id));
    }else{
      feature->vacate(siteId);
      feature->occupy(siteId);

      walker.setDwellTime(feature->getDwellTime(walker_id));
      walker.setPotentialSite(feature->pickNewSiteId(walker_id));
    }

    ++iteration_;
    if(iteration_ > iteration_threshold_){
      if(iteration_threshold_min_!=constants::inf_iterations){
        if(coarseGrain_(siteToHopToId)){
          iteration_threshold_ = iteration_threshold_min_;
        }else{
          iteration_threshold_*=2;
        }
      }
      iteration_ = 0;
    }
  }

  /****************************************************************************
   * Internal Private Functions
   ****************************************************************************/

  bool KMC_CoarseGrainSystem::coarseGrain_(int siteId){
    BasinExplorer basin_explorer;
    auto basin_site_ids = basin_explorer.findBasin(*sites_,*clusters_,siteId);

    double internal_time_limit = getInternalTimeLimit_(basin_site_ids);

    if( sitesSatisfyEquilibriumCondition_(basin_site_ids, internal_time_limit) ){
      auto sites_and_clusters = getClustersOfSites(basin_site_ids);
      auto number_clusters = countUniqueClusters(sites_and_clusters);

      if(number_clusters==1 &&
          sites_and_clusters.begin()->second==constants::unassignedId)
      {
        createCluster_(basin_site_ids,internal_time_limit);
        return true;
      }else if(number_clusters!=1){
        // Joint clusters and sites to an existing cluster
        int favored_clusterId = getFavoredClusterId(sites_and_clusters);
        mergeSitesAndClusters_(sites_and_clusters,favored_clusterId);
        return true;
      }
    }
    return false;
  }

  size_t countUniqueClusters(const unordered_map<int,int> & sites_and_clusters){
    set<int> clusters;
    for(auto site_and_cluster : sites_and_clusters){
      clusters.insert(site_and_cluster.second);
    }
    return clusters.size();
  }

  int getFavoredClusterId(unordered_map<int,int> sites_and_clusters){
    int clusterId = constants::unassignedId;
    for(auto site_and_cluster : sites_and_clusters){
      if(site_and_cluster.second != constants::unassignedId){
        if(clusterId==constants::unassignedId || 
            site_and_cluster.second < clusterId){
          clusterId = site_and_cluster.second;
        }
      }
    }
    return clusterId;
  }

  // The first int is the site id the second int is the cluster id 
  unordered_map<int,int> KMC_CoarseGrainSystem::getClustersOfSites(const vector<int> & siteIds){
    unordered_map<int,int> sites_and_clusters;
    for(auto siteId : siteIds){
      if(sites_->partOfCluster(siteId)){
        sites_and_clusters[siteId]= sites_->getClusterIdOfSite(siteId);
      }else{
        sites_and_clusters[siteId]= constants::unassignedId;
      }
    }
    return sites_and_clusters;
  }

  int KMC_CoarseGrainSystem::createCluster_(vector<int> siteIds, double internal_time_limit) {
    LOG("Creating cluster from vector of sites", 1);

    KMC_Cluster cluster;
    cluster.setConvergenceMethod(KMC_Cluster::Method::converge_by_tolerance);
    cluster.setConvergenceTolerance(0.001);
    vector<KMC_Site> sites;
    for (auto siteId : siteIds){
      sites.push_back(sites_->getKMC_Site(siteId));
    }
    cluster.addSites(sites);
    cluster.updateProbabilitiesAndTimeConstant();

    double cluster_time_const = cluster.getTimeConstant();
    // Cut the resolution in half from what it would otherwise be otherwise not worth doing
    double res = cluster_time_const/(2*internal_time_limit);
    double allowed_resolution = cluster_time_const/time_resolution_;
    double chosen_resolution = res;

    // The coarser the resolution is the better
    if(allowed_resolution <  chosen_resolution) chosen_resolution=allowed_resolution;

    if(chosen_resolution<2.0) chosen_resolution=2.0;
   
    cluster.setResolution(chosen_resolution);
    if (seed_set_) {
      cluster.setRandomSeed(seed_);
      ++seed_;
    }
    clusters_->addKMC_Cluster(cluster);

    for(auto siteId : siteIds){
      sites_->setClusterId(siteId,cluster.getId());  
      topology_features_[siteId] = &(clusters_->getKMC_Cluster(cluster.getId()));
    }

    auto sitesFoundInCluster = cluster.getSiteIdsInCluster();

    return cluster.getId();
  }

  void KMC_CoarseGrainSystem::mergeSitesAndClusters_( unordered_map<int,int> sites_and_clusters,int favoredClusterId) {

    LOG("Merging sites to cluster", 1);
    vector<KMC_Site> isolated_sites;
    unordered_set<int> cluster_ids;

    for (auto site_and_cluster : sites_and_clusters) { 
      if(site_and_cluster.second != favoredClusterId){ 
        if (site_and_cluster.second == constants::unassignedId) {
          isolated_sites.push_back(sites_->getKMC_Site(site_and_cluster.first));
        } else {
          cluster_ids.insert(site_and_cluster.second);
        }
        topology_features_[site_and_cluster.first] = &(clusters_->getKMC_Cluster(favoredClusterId));
        sites_->setClusterId(site_and_cluster.first,favoredClusterId);
      }
    }
    clusters_->getKMC_Cluster(favoredClusterId).addSites(isolated_sites);
    clusters_->getKMC_Cluster(favoredClusterId).updateProbabilitiesAndTimeConstant();
    for(auto clusterId : cluster_ids ){
      clusters_->getKMC_Cluster(favoredClusterId).migrateSitesFrom(clusters_->getKMC_Cluster(clusterId));
      clusters_->erase(clusterId);
    }

  }

  double KMC_CoarseGrainSystem::getExternalTimeLimit_(const vector<int> & siteIds ){
    LOG("Getting the external time limit of a cluster", 1);
    unordered_set<int> internal_sites;
    for(const int & site_id : siteIds){
      internal_sites.insert(site_id);
    }

    double max_rate_off = 0; 
    for(const int & site_id : siteIds){
      KMC_Site & site = sites_->getKMC_Site(site_id);
      const unordered_map<int,double *> neigh_and_rates = site.getNeighborsAndRatesConst();
      for( const pair<int,double *> & neigh_and_rate : neigh_and_rates){
        if(internal_sites.count(neigh_and_rate.first)==0){
          if(*(neigh_and_rate.second) > max_rate_off){
            max_rate_off = *(neigh_and_rate.second);
          }
        }
      }
    }
    return max_rate_off;
  }

double KMC_CoarseGrainSystem::getInternalTimeLimit_(vector<int> siteIds ){
  LOG("Getting the internal time limit of a cluster", 1);

  auto nodes = convertSitesToEmptySharedNodes(siteIds);

  unordered_map<int, weak_ptr<GraphNode<string>>> nodes_weak;
  for (auto node_iter : nodes) nodes_weak[node_iter.first] = node_iter.second;

  auto edges = convertSitesOutgoingRatesToTimeSharedWeightedEdges<vector<shared_ptr<Edge>>>(
      *sites_,
      siteIds);

  list<weak_ptr<Edge>> edges_weak(edges.begin(), edges.end());

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

// Its not worth creating a cluster unless the time is at least cut in half
// And it is not allowed if the sample time is smaller than than the simulated
// time of the cluster. The cluster has to be updated at a minimum once between
// each measurment (time_resolution). If this is not done the noise will not
// correctly show up in the data.  
// The number 25 is the ratio needed between hops within the cluster to hops
// outside of the cluster in order to see performance gains.
bool KMC_CoarseGrainSystem::sitesSatisfyEquilibriumCondition_(
    vector<int> siteIds, double maxtime) {

  LOG("Checking if sites satisfy equilibrium condition", 1);
  double timeConstant = getTimeConstantFromSitesToNeighbors_(siteIds);
  double time_to_traverse_cluster = maxtime*minimum_coarse_graining_resolution_;
  return timeConstant > time_to_traverse_cluster*performance_ratio_ && time_to_traverse_cluster< time_resolution_;// && ratio>25;
}

double KMC_CoarseGrainSystem::getTimeConstantFromSitesToNeighbors_(
   const vector<int> & siteIds) const {

  LOG("Get the minimum time constant", 1);
  set<int> internalSiteIds(siteIds.begin(), siteIds.end());

  double sumRates = 0.0;
  for (const int & siteId : siteIds) {
    vector<int> neighborSiteIds = sites_->getSiteIdsOfNeighbors(siteId);
    for (const int & neighId : neighborSiteIds) {
      if (!internalSiteIds.count(neighId)) {
        sumRates+= sites_->getRateToNeighborOfSite(siteId,neighId);
      }
    }
  }
  if (sumRates == 0.0){
    return 0.0;
  }
  return 1.0/sumRates;
}

unordered_map<int,vector<int>> KMC_CoarseGrainSystem::getClusters(){
  return clusters_->getSiteIdsOfClusters();
}

unordered_map<int,double> KMC_CoarseGrainSystem::getResolutionOfClusters(){
  return clusters_->getResolutionOfClusters();
}

unordered_map<int,double> KMC_CoarseGrainSystem::getTimeIncrementOfClusters(){
  return clusters_->getTimeIncrementOfClusters();
}

int KMC_CoarseGrainSystem::getFavoredClusterId_(vector<int> siteIds) {

  LOG("Getting the favored cluster Id", 1);
  int favoredClusterId = constants::unassignedId;
  for (auto siteId : siteIds) {
    int clusterId = sites_->getClusterIdOfSite(siteId);
    if (favoredClusterId == constants::unassignedId) {
      favoredClusterId = clusterId;
    } else if (clusterId != constants::unassignedId &&
               clusterId < favoredClusterId) {
      favoredClusterId = clusterId;
    }
  }
  return favoredClusterId;
}


}
