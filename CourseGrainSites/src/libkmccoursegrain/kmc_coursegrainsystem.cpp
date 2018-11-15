
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <chrono>

#include "../../include/kmccoursegrain/kmc_coursegrainsystem.hpp"
#include "../../include/kmccoursegrain/kmc_constants.hpp"
#include "../../include/kmccoursegrain/kmc_particle.hpp"

#include "topologyfeatures/kmc_topology_feature.hpp"
#include "topologyfeatures/kmc_cluster.hpp"
#include "topologyfeatures/kmc_site.hpp"
#include "log.hpp"
#include "kmc_basin_explorer.hpp"
#include "kmc_graph_library_adapter.hpp"

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

namespace kmccoursegrain {
/****************************************************************************
 * Private Internal Function Declarations
 ****************************************************************************/

bool compare(const pair<int,double> &x, const pair<int,double> &y){
  return x.second>y.second;
}

unordered_map<int, shared_ptr<GraphNode<string>>> createNode_(int siteIds);

size_t countUniqueClusters(unordered_map<int,int> sites_and_clusters);
int getFavoredClusterId(unordered_map<int,int> sites_and_clusters);

/****************************************************************************
 * Public Facing Functions
 ****************************************************************************/

void KMC_CourseGrainSystem::initializeSystem(unordered_map<int, unordered_map<int, double*>> ratesOfAllSites) {

  LOG("Initializeing system", 1);

  for (auto it = ratesOfAllSites.begin(); it != ratesOfAllSites.end(); ++it) {
    KMC_Site site;
    site.setId(it->first);
    site.setRatesToNeighbors(it->second);
    if (seed_set_) {
      site.setRandomSeed(seed_);
      ++seed_;
    }
    sites_.addKMC_Site(site);
    topology_features_[it->first] = &(sites_.getKMC_Site(it->first));
  }

}

void KMC_CourseGrainSystem::initializeParticles(vector<KMC_Particle>& particles) {

  LOG("Initializeing particles", 1);
  
  if (topology_features_.size() == 0) {
    throw runtime_error(
        "You must first initialize the system before you "
        "can initialize the particles");
  }

  for ( size_t index = 0; index<particles.size(); ++index){
//  for (KMC_Particle & particle : particles) {
    auto siteId = particles.at(index).getIdOfSiteCurrentlyOccupying();
    if (siteId == constants::unassignedId) {
      throw runtime_error(
          "You must first place the particle on a known site"
          " before the particle can be initialized.");
    }
    topology_features_[siteId]->occupy();

    int newId = topology_features_[siteId]->pickNewSiteId();
    auto hopTime = topology_features_[siteId]->getDwellTime();
    particles.at(index).setDwellTime(hopTime);
    particles.at(index).setPotentialSite(newId);
  }
}

void KMC_CourseGrainSystem::setMinCourseGrainIterationThreshold(int threshold_min) {
  LOG("Setting minimum course graining threshold", 1);
  iteration_threshold_min_ = threshold_min;
  iteration_threshold_ = threshold_min;
}

void KMC_CourseGrainSystem::setRandomSeed(const unsigned long seed) {
  if (topology_features_.size() != 0) {
    throw runtime_error(
        "For the random seed to have an affect, it must be "
        "set before initializeSystem is called");
  }
  seed_ = seed;
  seed_set_ = true;
}

void KMC_CourseGrainSystem::removeParticleFromSystem(KMC_Particle& particle) {
  LOG("Particle is being removed from system", 1);
  auto siteId = particle.getIdOfSiteCurrentlyOccupying();
  sites_.vacate(siteId);
}

int KMC_CourseGrainSystem::getClusterIdOfSite(int siteId) {
  return sites_.getClusterIdOfSite(siteId);
}

void KMC_CourseGrainSystem::hop(KMC_Particle & particle) {
  LOG("Particle is hopping in system", 1);
  auto siteId = particle.getIdOfSiteCurrentlyOccupying();
  int siteToHopToId = particle.getPotentialSite();

  KMC_TopologyFeature * feature = topology_features_[siteId];
  KMC_TopologyFeature * feature_to_hop_to = topology_features_[siteToHopToId];

  int newId;
  double hopTime;
  
  if(!feature_to_hop_to->isOccupied(siteToHopToId)){
    LOG("Hopping to " + to_string(siteToHopToId) + " site", 1);
    feature->vacate(siteId);
    feature_to_hop_to->occupy(siteToHopToId);

    newId   = feature_to_hop_to->pickNewSiteId();
    hopTime = feature_to_hop_to->getDwellTime();
 
    particle.occupySite(siteToHopToId);
    particle.setDwellTime(hopTime);
    particle.setPotentialSite(newId);
  }else{
//    siteToHopToId = siteId;
//    feature->vacate(siteId);
//    feature_to_hop_to->occupy(siteId);

    newId   = feature_to_hop_to->pickNewSiteId();
    hopTime = feature_to_hop_to->getDwellTime();
    
    particle.setDwellTime(hopTime);
    particle.setPotentialSite(newId);
  }

  ++iteration_;
  if(iteration_ > iteration_threshold_){
    if(courseGrain_(siteToHopToId)){
      iteration_threshold_ = iteration_threshold_min_;
    }else{
      iteration_threshold_*=2;
    }
  }
}

/****************************************************************************
 * Internal Private Functions
 ****************************************************************************/

bool KMC_CourseGrainSystem::courseGrain_(int siteId){
  BasinExplorer basin_explorer;
  auto basin_site_ids = basin_explorer.findBasin(sites_,siteId);

  double internal_time_limit = getInternalTimeLimit_(basin_site_ids);
  if( sitesSatisfyEquilibriumCondition_(basin_site_ids, internal_time_limit) ){

    auto sites_and_clusters = getClustersOfSites(basin_site_ids);
    auto number_clusters = countUniqueClusters(sites_and_clusters);

    if(number_clusters==1 &&
       sites_and_clusters.begin()->second==constants::unassignedId)
    {
      cout << "Creating cluster " << endl;
      createCluster_(basin_site_ids,internal_time_limit);
      return true;
    }else if(number_clusters!=1){
      cout << "Merging cluster " << endl;
      // Joint clusters and sites to an existing cluster
      int favored_clusterId = getFavoredClusterId(sites_and_clusters);
      mergeSitesAndClusters_(sites_and_clusters,favored_clusterId);
      return true;
    }
  }
  return false;
}

size_t countUniqueClusters(unordered_map<int,int> sites_and_clusters){
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
unordered_map<int,int> KMC_CourseGrainSystem::getClustersOfSites(vector<int> siteIds){
  unordered_map<int,int> sites_and_clusters;
  for(auto siteId : siteIds){
    if(sites_.partOfCluster(siteId)){
      sites_and_clusters[siteId]= sites_.getClusterIdOfSite(siteId);
    }else{
      sites_and_clusters[siteId]= constants::unassignedId;
    }
  }
  return sites_and_clusters;
}

int KMC_CourseGrainSystem::createCluster_(vector<int> siteIds, double internal_time_limit) {
  LOG("Creating cluster from vector of sites", 1);

  KMC_Cluster cluster;
  cluster.setConvergenceMethod(KMC_Cluster::Method::converge_by_tolerance);
  cluster.setConvergenceTolerance(0.001);
  vector<KMC_Site> sites;
  for (auto siteId : siteIds){
    sites.push_back(sites_.getKMC_Site(siteId));
  }
  cluster.addSites(sites);

  double cluster_time_const = cluster.getTimeConstant();
  // Cut the resolution in half from what it would otherwise be otherwise not worth doing
  int res = static_cast<int>(floor(cluster_time_const/(2*internal_time_limit)));
  int chosen_resolution = res;
  if(res==0) chosen_resolution=1;
  if(clusterResolution_<res) chosen_resolution = clusterResolution_; 
  cluster.setResolution(chosen_resolution);
  if (seed_set_) {
    cluster.setRandomSeed(seed_);
    ++seed_;
  }
  clusters_[cluster.getId()] = cluster;

  for(auto siteId : siteIds){
    sites_.setClusterId(siteId,cluster.getId());  
    topology_features_[siteId] = &(clusters_[cluster.getId()]);
  }

  auto sitesFoundInCluster = cluster.getSiteIdsInCluster();

  return cluster.getId();
}

void KMC_CourseGrainSystem::mergeSitesAndClusters_( unordered_map<int,int> sites_and_clusters,int favoredClusterId) {

  LOG("Merging sites to cluster", 1);
  vector<KMC_Site> isolated_sites;
  unordered_set<int> cluster_ids;

  for (auto site_and_cluster : sites_and_clusters) { 

    if(site_and_cluster.second != favoredClusterId){ 
      if (site_and_cluster.second == constants::unassignedId) {
        isolated_sites.push_back(sites_.getKMC_Site(site_and_cluster.first));
      } else {
        cluster_ids.insert(site_and_cluster.second);
      }
      topology_features_[site_and_cluster.first] = &(clusters_[favoredClusterId]);
    }
  }
  clusters_[favoredClusterId].addSites(isolated_sites);
  for(auto clusterId : cluster_ids ){
    clusters_[favoredClusterId].migrateSitesFrom(clusters_[clusterId]);
    clusters_.erase(clusters_.find(clusterId));
  }

}

double KMC_CourseGrainSystem::getInternalTimeLimit_(vector<int> siteIds ){
  LOG("Getting the internal time limit of a cluster", 1);

  auto nodes = convertSitesToEmptySharedNodes(siteIds);

  unordered_map<int, weak_ptr<GraphNode<string>>> nodes_weak;
  for (auto node_iter : nodes) nodes_weak[node_iter.first] = node_iter.second;

  auto edges = convertSitesOutgoingRatesToTimeSharedWeightedEdges<vector<shared_ptr<Edge>>>(
      sites_,
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

// Its not worth creating a cluster unless the time being at least cut in half
bool KMC_CourseGrainSystem::sitesSatisfyEquilibriumCondition_(
    vector<int> siteIds, double maxtime) {

  LOG("Checking if sites satisfy equilibrium condition", 1);
  auto minTimeConstant = getMinimumTimeConstantFromSitesToNeighbors_(siteIds);
  return minTimeConstant > (maxtime*minimum_course_graining_resolution_);
}

double KMC_CourseGrainSystem::getMinimumTimeConstantFromSitesToNeighbors_(
    vector<int> siteIds) {

  LOG("Get the minimum time constant", 1);
  set<int> internalSiteIds(siteIds.begin(), siteIds.end());

  double minimumTimeConstant = -1.0;
  bool initialized = false;
  for (auto siteId : siteIds) {
    auto neighborSiteIds = sites_.getSiteIdsOfNeighbors(siteId);
    for (auto neighId : neighborSiteIds) {
      if (!internalSiteIds.count(neighId)) {
        auto timeConstant = 1.0 / sites_.getRateToNeighborOfSite(siteId,neighId);
        if (!initialized) {
          initialized = true;
          minimumTimeConstant = timeConstant;
        } else if (timeConstant < minimumTimeConstant) {
          minimumTimeConstant = timeConstant;
        }
      }
    }
  }
  return minimumTimeConstant;
}

vector<vector<int>> KMC_CourseGrainSystem::getClusters(){
  vector<vector<int>> clusters;
  for(auto cluster : clusters_){
    clusters.push_back(cluster.second.getSiteIdsInCluster());
  }
  return clusters;
}

int KMC_CourseGrainSystem::getFavoredClusterId_(vector<int> siteIds) {

  LOG("Getting the favored cluster Id", 1);
  int favoredClusterId = constants::unassignedId;
  for (auto siteId : siteIds) {
    int clusterId = sites_.getClusterIdOfSite(siteId);
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
