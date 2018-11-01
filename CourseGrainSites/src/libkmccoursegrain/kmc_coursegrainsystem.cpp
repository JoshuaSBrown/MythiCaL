#include "../../include/kmccoursegrain/kmc_coursegrainsystem.hpp"
#include "../../include/kmccoursegrain/kmc_constants.hpp"
#include "../../include/kmccoursegrain/kmc_particle.hpp"
#include "kmc_cluster.hpp"
#include "kmc_site.hpp"
#include "log.hpp"
#include <iostream>
#include <map>
#include <stdexcept>
#include <algorithm>

#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph.hpp"
#include "../../../UGLY/include/ugly/graph_algorithms.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"

using namespace std;
using namespace ugly;
using namespace ugly::graphalgorithms;

namespace kmccoursegrain {
/****************************************************************************
 * Private Internal Function Declarations
 ****************************************************************************/

/**
 * \brief Create a list of edges that will define a graph
 *
 * This function is designed to take a list of site ids and calcualte the
 * edges that connect them.
 *
 * \param[in] siteIds - a vector of the siteIds that will potentially make a
 * cluster
 * \param[in] sites - map of pointers to the sites
 *
 * \return list of edges defining the connections between the sites
 **/
list<shared_ptr<Edge>> createEdges_(unordered_map<int, SitePtr> sites_,
                                    vector<int> siteIds);
list<shared_ptr<Edge>> createEdges_(unordered_map<int, SitePtr> sites_,
                                    unordered_map<int,int> siteIds);


bool compare(const pair<int,int> &x, const pair<int,int> &y){
  return x.second>y.second;
}
/**
 * \brief Creates a map of GraphNodes that define the vertices of a graph
 *
 * The site ids that will potentially make a cluster are passed in and turned
 * into vertices.
 *
 * \param[in] siteIds vector of the ids that will make up the graph vertices
 *
 * \return a map of nodes where the first int is the vertex id and the second
 * consists of information stored in the node or vertex in the form of a
 * graphnode
 **/
map<int, shared_ptr<GraphNode<string>>> createNodes_(vector<int> siteIds);
map<int, shared_ptr<GraphNode<string>>> createNodes_(unordered_map<int,int> siteIds);

/****************************************************************************
 * Public Facing Functions
 ****************************************************************************/

void KMC_CourseGrainSystem::initializeSystem(
    map<int const, map<int const, double*>> ratesOfAllSites) {

  LOG("Initializeing system", 1);

  for (auto it = ratesOfAllSites.begin(); it != ratesOfAllSites.end(); ++it) {
    auto site = shared_ptr<KMC_Site>(new KMC_Site);
    site->setId(it->first);
    site->setRatesToNeighbors(it->second);
    if (seed_set_) {
      site->setRandomSeed(seed_);
      ++seed_;
    }
    sites_[it->first] = site;
  }
}

void KMC_CourseGrainSystem::initializeParticles(vector<ParticlePtr> particles) {

  LOG("Initializeing particles", 1);
  
  if (sites_.size() == 0) {
    throw runtime_error(
        "You must first initialize the system before you "
        "can initialize the particles");
  }

  for (auto particle : particles) {
    auto siteId = particle->getIdOfSiteCurrentlyOccupying();
    if (siteId == constants::unassignedId) {
      throw runtime_error(
          "You must first place the particle on a known site"
          " before the particle can be initialized.");
    }
    sites_[siteId]->occupySite();

    int newId = sites_[siteId]->pickNewSiteId();
    auto hopTime = sites_[siteId]->getDwellTime();
    particle->setDwellTime(hopTime);
    particle->setPotentialSite(newId);
  }
}

void KMC_CourseGrainSystem::setCourseGrainIterationThreshold(int threshold) {
  LOG("Setting threshold", 1);
  iteration_threshold_ = threshold;
}

void KMC_CourseGrainSystem::setRandomSeed(const unsigned long seed) {
  if (sites_.size() != 0) {
    throw runtime_error(
        "For the random seed to have an affect, it must be "
        "set before initializeSystem is called");
  }
  seed_ = seed;
  seed_set_ = true;
}

void KMC_CourseGrainSystem::removeParticleFromSystem(ParticlePtr particle) {
  LOG("Particle is being removed from system", 1);
  auto siteId = particle->getIdOfSiteCurrentlyOccupying();
  sites_[siteId]->vacateSite();
}

int KMC_CourseGrainSystem::getClusterIdOfSite(int siteId) {
  return sites_[siteId]->getClusterId();
}

void KMC_CourseGrainSystem::hop(ParticlePtr particle) {
  LOG("Particle is hopping in system", 1);
  auto siteId = particle->getIdOfSiteCurrentlyOccupying();
  int siteToHopTo = particle->getPotentialSite();

  ++iteration_;

  if (sites_[siteToHopTo]->siteIsOccupied()) {
    LOG("Site " + to_string(siteToHopTo) + " is occupied", 1);
    if (sites_[siteId]->partOfCluster()) {
      
      auto clusterId = sites_[siteId]->getClusterId();
      int newId = clusters_[clusterId]->pickNewSiteId();
      auto hopTime = clusters_[clusterId]->getDwellTime();
      particle->setDwellTime(hopTime);
      particle->setPotentialSite(newId);

    } else {

      int newId = sites_[siteId]->pickNewSiteId();
      auto hopTime = sites_[siteId]->getDwellTime();

      particle->setDwellTime(hopTime);
      particle->setPotentialSite(newId);
    }

  } else {

    LOG("Hopping to " + to_string(siteToHopTo) + " site", 1);
    if (sites_[siteToHopTo]->partOfCluster()) {
      auto clusterId = sites_[siteToHopTo]->getClusterId();
      sites_[siteId]->vacateSite();
      sites_[siteToHopTo]->occupySite();

      int newId = clusters_[clusterId]->pickNewSiteId();
      auto hopTime = clusters_[clusterId]->getDwellTime();
      particle->occupySite(siteToHopTo);
      particle->setDwellTime(hopTime);
      particle->setPotentialSite(newId);
      //courseGrainSiteIfNeeded_(particle);

    } else {
      sites_[siteId]->vacateSite();
      sites_[siteToHopTo]->occupySite();

      int newId = sites_[siteToHopTo]->pickNewSiteId();
      auto hopTime = sites_[siteToHopTo]->getDwellTime();
      particle->occupySite(siteToHopTo);
      particle->setDwellTime(hopTime);
      particle->setPotentialSite(newId);

      sites_visited_.insert(newId);
      //courseGrainSiteIfNeeded_(particle);
    }
    if(sites_[siteToHopTo]->getVisitFrequency()>max_sample_frequency_) {
      max_sample_frequency_=sites_[siteId]->getVisitFrequency();
    }
  }

  if(iteration_ > iteration_threshold_){
    auto relevantSites = filterSites_();
    breakIntoIslands_(relevantSites);
    max_sample_frequency_ = max_frequency_;
    iteration_ = 0;
  }
}

/****************************************************************************
 * Internal Private Functions
 ****************************************************************************/

vector<vector<int>> KMC_CourseGrainSystem::breakIntoIslands_(unordered_map<int,int> relevantSites){

  auto edges = createEdges_(sites_,relevantSites);
  auto nodes = createNodes_(relevantSites);

//  cout << "number of edges " << edges.size() << endl;
  list<weak_ptr<Edge>> edges_weak(edges.begin(), edges.end());
  map<int, weak_ptr<GraphNode<string>>> nodes_weak;
  for (auto map_iter : nodes) nodes_weak[map_iter.first] = map_iter.second;

  auto graph_ptr = shared_ptr<Graph<string>>(new Graph<string>(edges_weak, nodes_weak));
  auto initial_islands = findSubGraphs(*graph_ptr); 
  vector<vector<int>> islands;

//  cout << "Number of relevant sites " << relevantSites.size() << endl;
//  cout << "Number of islands " << initial_islands.size() << endl;
  for ( auto island : initial_islands ){
//    cout << "island " << endl;
  // Sort the islands from highest frequency to lowest
    vector<pair<int,int>> island_frequencies;
    for( auto v : island){
      island_frequencies.push_back(pair<int,int>(v,relevantSites[v]));
    }
    sort(island_frequencies.begin(),island_frequencies.end(),compare);

//    cout << "Island " << endl;
//    for(auto pr : island_frequencies){
//      cout << pr.first << " " << pr.second << endl;
//    }
    if(island_frequencies.size()>5){
      for( auto v : island ){
        sampled_sites_.insert(v);
      }
    }else{
      for ( size_t count =1 ;count< island.size();++count){

        vector<int> tempIsland;
        for( size_t ind = 0; ind<=count;++ind){
          tempIsland.push_back(island_frequencies.at(ind).first);
        }

        auto edges2 = createEdges_(sites_,tempIsland);
        auto nodes2 = createNodes_(tempIsland);

        list<weak_ptr<Edge>> edges_weak2(edges2.begin(), edges2.end());
        map<int, weak_ptr<GraphNode<string>>> nodes_weak2;
        for (auto map_iter : nodes2) nodes_weak2[map_iter.first] = map_iter.second;

        auto graph_ptr2 = shared_ptr<Graph<string>>(new Graph<string>(edges_weak2, nodes_weak2));
        if(isSingleNetwork(*graph_ptr2)){
          double internal_time_limit = getInternalTimeLimit_(tempIsland);
          if( sitesSatisfyEquilibriumCondition_(tempIsland, internal_time_limit) ){
            createCluster_(tempIsland,internal_time_limit);
            // Remove sites from consideration

            for( auto v : island ){
              sampled_sites_.insert(v);
            }
            
            break;
          }
        }

      }

    } 
  }

  return islands;

}

unordered_map<int,int> KMC_CourseGrainSystem::filterSites_(){
  unordered_map<int,int> high_frequency_sites;

  unordered_set<int> important_sites;
  for( auto siteId : sites_visited_ ){
    auto freq = sites_[siteId]->getVisitFrequency();
    if(freq > max_sample_frequency_/ratio_threshold_relevant_sites_){
      if(sampled_sites_.count(siteId)==0){
        important_sites.insert(siteId);
        high_frequency_sites[siteId]= freq;
      }
    } 
  }
  return high_frequency_sites;
}

void KMC_CourseGrainSystem::createCluster_(vector<int> siteIds, double internal_time_limit) {
  LOG("Creating cluster from vector of sites", 1);
  auto cluster_ptr = shared_ptr<KMC_Cluster>(new KMC_Cluster());
  vector<SitePtr> sites;
  for (auto siteId : siteIds) sites.push_back(sites_[siteId]);

  cluster_ptr->addSites(sites);

  double cluster_dwell = cluster_ptr->getDwellTime();
  int res = static_cast<int>(floor(cluster_dwell/internal_time_limit));
  int chosen_resolution = res;
  if(res==0) chosen_resolution=1;
  if(clusterResolution_<res) chosen_resolution = clusterResolution_; 

  cout << "Chosen resolution " << chosen_resolution << endl;
  cluster_ptr->setResolution(chosen_resolution);
  if (seed_set_) {
    cluster_ptr->setRandomSeed(seed_);
    ++seed_;
  }
  clusters_[cluster_ptr->getId()] = cluster_ptr;
}

void KMC_CourseGrainSystem::mergeSitesToCluster_(vector<int> siteIds,
                                                 int favoredClusterId) {

  LOG("Merging sites to cluster", 1);
  vector<SitePtr> isolated_sites;
  for (auto siteId : siteIds) {
    int clusterId = sites_[siteId]->getClusterId();
    if (clusterId == constants::unassignedId) {
      isolated_sites.push_back(sites_[siteId]);
    } else if (clusterId != favoredClusterId) {
      clusters_[favoredClusterId]->migrateSitesFrom(clusters_[clusterId]);
      clusters_.erase(clusterId);
    }
  }
  clusters_[favoredClusterId]->addSites(isolated_sites);
}

double KMC_CourseGrainSystem::getInternalTimeLimit_(vector<int> siteIds){
  LOG("Getting the internal time limit of a cluster", 1);
  auto edges = createEdges_(sites_, siteIds);
  auto nodes = createNodes_(siteIds);
  list<weak_ptr<Edge>> edges_weak(edges.begin(), edges.end());
  map<int, weak_ptr<GraphNode<string>>> nodes_weak;
  for (auto map_iter : nodes) nodes_weak[map_iter.first] = map_iter.second;

  auto graph_ptr =
      shared_ptr<Graph<string>>(new Graph<string>(edges_weak, nodes_weak));

  map<pair<int, int>, double> verticesAndtimes =
      maxMinimumDistanceBetweenEveryVertex<string>(*graph_ptr);

  double maxtime = 0.0;
  for (auto verticesAndTime : verticesAndtimes) {
    if (verticesAndTime.second > maxtime) maxtime = verticesAndTime.second;
  }
  return maxtime;
}

bool KMC_CourseGrainSystem::sitesSatisfyEquilibriumCondition_(
    vector<int> siteIds, double maxtime) {

  LOG("Checking if sites satisfy equilibrium condition", 1);
  auto minTimeConstant = getMinimumTimeConstantFromSitesToNeighbors_(siteIds);
//  cout << "Min time " << minTimeConstant << " maxtime " << maxtime << endl;
  return minTimeConstant > (maxtime);
}

double KMC_CourseGrainSystem::getMinimumTimeConstantFromSitesToNeighbors_(
    vector<int> siteIds) {

  LOG("Get the minimum time constant", 1);
  set<int> internalSiteIds(siteIds.begin(), siteIds.end());

  double minimumTimeConstant = -1.0;
  bool initialized = false;
  for (auto siteId : siteIds) {
    auto neighborSiteIds = sites_[siteId]->getNeighborSiteIds();
    for (auto neighId : neighborSiteIds) {
      if (!internalSiteIds.count(neighId)) {
        auto timeConstant = 1.0 / sites_[siteId]->getRateToNeighbor(neighId);
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

map<int, shared_ptr<GraphNode<string>>> createNodes_(vector<int> siteIds) {
  LOG("Creating nodes", 1);
  map<int, shared_ptr<GraphNode<string>>> nds;
  for (auto siteId : siteIds) {
    if (nds.count(siteId) == 0) {
      auto gn = shared_ptr<GraphNode<string>>(new GraphNode<string>(""));
      nds[siteId] = gn;
    }
  }
  return nds;
}

map<int, shared_ptr<GraphNode<string>>> createNodes_(unordered_map<int,int> siteIds) {
  LOG("Creating nodes", 1);
  map<int, shared_ptr<GraphNode<string>>> nds;
  for (auto siteId : siteIds) {
    if (nds.count(siteId.first) == 0) {
      auto gn = shared_ptr<GraphNode<string>>(new GraphNode<string>(""));
      nds[siteId.first] = gn;
    }
  }
  return nds;
}

list<shared_ptr<Edge>> createEdges_(map<int, SitePtr> sites_,
                                    vector<int> siteIds) {

  LOG("Creating graph", 1);
  // Need edges to be directed and weighted
  list<shared_ptr<Edge>> edges;
  for (auto siteId : siteIds) {
    for (auto potential_neighId : siteIds) {
      if (sites_[siteId]->isNeighbor(potential_neighId)) {
        auto time = 1.0 / sites_[siteId]->getRateToNeighbor(potential_neighId);
        auto edge = shared_ptr<EdgeDirectedWeighted>(
            new EdgeDirectedWeighted(siteId, potential_neighId, time));

        edges.push_back(edge);
      }
    }
  }

  return edges;
}

list<shared_ptr<Edge>> createEdges_(map<int, SitePtr> sites_,
                                    unordered_map<int,int> siteIds) {

  LOG("Creating graph", 1);
  // Need edges to be directed and weighted
  list<shared_ptr<Edge>> edges;
  for (auto siteId : siteIds) {
    for (auto potential_neighId : siteIds) {
      if (sites_[siteId.first]->isNeighbor(potential_neighId.first)) {
        auto time = 1.0 / sites_[siteId.first]->getRateToNeighbor(potential_neighId.first);
        auto edge = shared_ptr<EdgeDirectedWeighted>(
            new EdgeDirectedWeighted(siteId.first, potential_neighId.first, time));

        edges.push_back(edge);
      }
    }
  }

  return edges;
}

int KMC_CourseGrainSystem::getFavoredClusterId_(vector<int> siteIds) {

  LOG("Getting the favored cluster Id", 1);
  int favoredClusterId = constants::unassignedId;
  for (auto siteId : siteIds) {
    int clusterId = sites_[siteId]->getClusterId();
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
