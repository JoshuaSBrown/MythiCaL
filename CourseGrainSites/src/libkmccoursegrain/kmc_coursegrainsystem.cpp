#include <stdexcept>
#include <iostream>
#include <map>
#include "log.hpp"
#include "kmc_cluster.hpp"
#include "kmc_site.hpp"
#include "../../include/kmccoursegrain/kmc_coursegrainsystem.hpp"
#include "../../include/kmccoursegrain/kmc_particle.hpp"
#include "../../include/kmccoursegrain/kmc_constants.hpp"

#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"
#include "../../../UGLY/include/ugly/graph_algorithms.hpp"

using namespace std;
using namespace ugly;
using namespace ugly::graphalgorithms;

namespace kmccoursegrain {
  /****************************************************************************
   * Private Internal Function Declarations
   ****************************************************************************/

  list<shared_ptr<Edge>> createEdges_(map<int,SitePtr> sites_, vector<int> siteIds);

  map<int,shared_ptr<GraphNode<string>>> createNodes_(vector<int> siteIds);

  /****************************************************************************
   * Public Facing Functions
   ****************************************************************************/

  void KMC_CourseGrainSystem::initializeSystem(
      map<int const,map<int const,double* >> ratesOfAllSites){

    LOG("Initializeing system",1);

    for(auto it = ratesOfAllSites.begin();it!=ratesOfAllSites.end();++it){
      auto site = shared_ptr<KMC_Site>(new KMC_Site);
      site->setId(it->first);
      site->setRatesToNeighbors(it->second);
      if(seed_set_) {
        site->setRandomSeed(seed_);
        ++seed_;
      }
      sites_[it->first] = move(site);
    }
  }

  void 
  KMC_CourseGrainSystem::initializeParticles(vector<ParticlePtr> particles){

    LOG("Initializeing particles",1);

    if(sites_.size()==0){
      throw runtime_error("You must first initialize the system before you "
          "can initialize the particles");
    }

    for(auto particle : particles){
      auto siteId = particle->getIdOfSiteCurrentlyOccupying();
      if(siteId == constants::unassignedId){
        throw runtime_error("You must first place the particle on a known site"
            " before the particle can be initialized.");
      }
      sites_[siteId]->occupySite();

      int newId = sites_[siteId]->pickNewSiteId();
      auto hopTime = sites_[siteId]->getDwellTime();

      particle->setDwellTime(hopTime);
      particle->setPotentialSite(newId);
    }
  }

  void KMC_CourseGrainSystem::setCourseGrainThreshold(int threshold){
    LOG("Setting threshold",1);
    courseGrainingThreshold_ = threshold;
  }

  void KMC_CourseGrainSystem::setRandomSeed(const unsigned long seed){
    if(sites_.size()!=0){
      throw runtime_error("For the random seed to have an affect, it must be "
          "set before initializeSystem is called");
    }
    seed_=seed;
    seed_set_=true;
  }

  void KMC_CourseGrainSystem::removeParticleFromSystem(ParticlePtr particle){
    LOG("Particle is being removed from system",1);
    auto siteId = particle->getIdOfSiteCurrentlyOccupying();
    sites_[siteId]->vacateSite();
  }

  int KMC_CourseGrainSystem::getClusterIdOfSite(int siteId){
    return sites_[siteId]->getClusterId();
  }

  void KMC_CourseGrainSystem::hop(ParticlePtr particle){

    LOG("Particle is hopping in system",1);
    auto siteId = particle->getIdOfSiteCurrentlyOccupying();
    int siteToHopTo = particle->getPotentialSite();

    if( sites_[siteToHopTo]->siteIsOccupied() ){
      LOG("Site "+to_string(siteToHopTo)+" is occupied",1);
      if( sites_[siteId]->partOfCluster() ){

        auto clusterId = sites_[siteId]->getClusterId();
        int newId = clusters_[clusterId]->pickNewSiteId();
        auto hopTime = clusters_[clusterId]->getDwellTime();
        particle->setDwellTime(hopTime);
        particle->setPotentialSite(newId);

      }else{

        int newId = sites_[siteId]->pickNewSiteId();
        auto hopTime = sites_[siteId]->getDwellTime();

        particle->setDwellTime(hopTime);
        particle->setPotentialSite(newId);
      }

    }else{

      if( sites_[siteToHopTo]->partOfCluster() ){
        LOG("Hopping to "+to_string(siteToHopTo)+" site",1);
        auto clusterId = sites_[siteToHopTo]->getClusterId();
        sites_[siteId]->vacateSite();
        sites_[siteToHopTo]->occupySite();

        int newId = clusters_[clusterId]->pickNewSiteId();
        auto hopTime = clusters_[clusterId]->getDwellTime();

        particle->occupySite(siteToHopTo,clusterId);
        particle->setDwellTime(hopTime);
        particle->setPotentialSite(newId);
        courseGrainSiteIfNeeded_(particle);

      }else{
        sites_[siteId]->vacateSite();
        sites_[siteToHopTo]->occupySite();

        int newId = sites_[siteToHopTo]->pickNewSiteId();
        auto hopTime = sites_[siteToHopTo]->getDwellTime();
        particle->occupySite(siteToHopTo,constants::unassignedId);
        particle->setDwellTime(hopTime);
        particle->setPotentialSite(newId);

        courseGrainSiteIfNeeded_(particle);

      }

    }

  }

  /****************************************************************************
   * Internal Private Functions
   ****************************************************************************/

  void 
  KMC_CourseGrainSystem::createCluster_(vector<int> siteIds){
    LOG("Creating cluster from vector of sites",1);
    auto cluster_ptr = shared_ptr<KMC_Cluster>(new KMC_Cluster());
    vector<SitePtr> sites;
    for(auto siteId : siteIds ) sites.push_back(sites_[siteId]);

    cluster_ptr->addSites(sites);

    if(seed_set_){
      cluster_ptr->setRandomSeed(seed_);
      ++seed_;
    }
    clusters_[cluster_ptr->getId()] = cluster_ptr;
  }

  void 
    KMC_CourseGrainSystem::mergeSitesToCluster_(
        vector<int> siteIds,int favoredClusterId){
  
    LOG("Merging sites to cluster",1);
    vector<SitePtr> isolated_sites;
    for(auto siteId : siteIds) {
      int clusterId = sites_[siteId]->getClusterId();
      if(clusterId==constants::unassignedId){
        isolated_sites.push_back(sites_[siteId]);
      }else if(clusterId!=favoredClusterId){
        clusters_[favoredClusterId]->migrateSitesFrom(clusters_[clusterId]);
        clusters_.erase(clusterId);
      }
    }
    clusters_[favoredClusterId]->addSites(isolated_sites);

  }

  bool 
    KMC_CourseGrainSystem::sitesSatisfyEquilibriumCondition_(
        vector<int> siteIds){

    LOG("Checking if sites satisfy equilibrium condition",1);
    auto edges = createEdges_(sites_,siteIds);
    auto nodes = createNodes_(siteIds);
    list<weak_ptr<Edge>> edges_weak(edges.begin(),edges.end());
    map<int,weak_ptr<GraphNode<string>>> nodes_weak;
    for(auto map_iter : nodes ) nodes_weak[map_iter.first] = map_iter.second;
 
    auto graph_ptr = shared_ptr<Graph<string>>\
                     (new Graph<string>(edges_weak,nodes_weak));

    map<pair<int,int>,double> verticesAndtimes = \
      maxMinimumDistanceBetweenEveryVertex<string>(*graph_ptr);

    double maxtime = 0.0;
    for( auto verticesAndTime : verticesAndtimes){
      if(verticesAndTime.second>maxtime) maxtime = verticesAndTime.second;
    }

    auto minTimeConstant = getMinimumTimeConstantFromSitesToNeighbors_(siteIds);
    return minTimeConstant>(2*maxtime);
  }
  
  double 
    KMC_CourseGrainSystem::getMinimumTimeConstantFromSitesToNeighbors_(
      vector<int> siteIds){
    
    LOG("Get the minimum time constant",1);
    set<int> internalSiteIds(siteIds.begin(),siteIds.end());

    double minimumTimeConstant=-1.0;
    bool initialized =false;
    for(auto siteId : siteIds){
      auto neighborSiteIds = sites_[siteId]->getNeighborSiteIds();  
      for( auto neighId : neighborSiteIds){
        if(!internalSiteIds.count(neighId)){
          auto timeConstant = 1.0/sites_[siteId]->getRateToNeighbor(neighId);
          if(!initialized){
            initialized=true;
            minimumTimeConstant = timeConstant;
          }else if(timeConstant<minimumTimeConstant){
            minimumTimeConstant = timeConstant;
          }
        }
      }
    }
    return minimumTimeConstant;
  }

  map<int,shared_ptr<GraphNode<string>>> createNodes_(vector<int> siteIds){
    LOG("Creating nodes",1);
    map<int,shared_ptr<GraphNode<string>>> nds;
    for(auto siteId : siteIds){
      if(nds.count(siteId)==0){
        auto gn = shared_ptr<GraphNode<string>>(new GraphNode<string>(""));
        nds[siteId] = gn;
      }
    }
    return nds;
  }

  list<shared_ptr<Edge>> createEdges_(map<int,SitePtr> sites_, vector<int> siteIds){

    LOG("Creating graph",1);
    // Need edges to be directed and weighted
    list<shared_ptr<Edge>> edges;
    for(auto siteId : siteIds){
      for(auto potential_neighId : siteIds){
        if(sites_[siteId]->isNeighbor(potential_neighId)){
          auto time = 1.0/sites_[siteId]->getRateToNeighbor(potential_neighId);
          auto edge = shared_ptr<EdgeDirectedWeighted>(new EdgeDirectedWeighted(siteId,potential_neighId,time));

          edges.push_back(edge);
        }
      }
    }

    return edges;
  }

  int KMC_CourseGrainSystem::getFavoredClusterId_(vector<int> siteIds){

    LOG("Getting the favored cluster Id",1);
    int favoredClusterId = constants::unassignedId;
    for( auto siteId : siteIds){
      int clusterId = sites_[siteId]->getClusterId();
      if(favoredClusterId==constants::unassignedId){
        favoredClusterId = clusterId;
      }else if( clusterId!=constants::unassignedId && clusterId<favoredClusterId){
        favoredClusterId = clusterId;
      }
    }
    return favoredClusterId;
  }

  void KMC_CourseGrainSystem::courseGrainSiteIfNeeded_(ParticlePtr& particle){

    LOG("Course graining sites if needed",1);

    auto memories = particle->getMemory();
 
    // Determine if the particle has made enough jumps to have a memory
    if(memories.size()>1){

      // Determine how many sites are above the threshold, stores the ones that
      // are over the threshold and appear in consecutive order from the most 
      // recently visited
      vector<int> relevantSites;
      for(auto memory : memories ){
        if(memory.at(2)>courseGrainingThreshold_){
          relevantSites.push_back(memory.at(0));
        }else{
          break;
        }
      }

      if(relevantSites.size()>1){

        bool satisfy = sitesSatisfyEquilibriumCondition_(relevantSites); 
        int favoredClusterId = getFavoredClusterId_(relevantSites); 
        if(satisfy){

          if(favoredClusterId==constants::unassignedId){
            // None of the sites are part of a cluster 
            createCluster_(relevantSites);
            int clusterId = sites_[relevantSites.at(0)]->getClusterId();
            // Remove the particles memory of visiting the second site
            for(int index = 1;index < static_cast<int>(relevantSites.size());++index){
              particle->removeMemory(relevantSites.at(index));
            }
            particle->resetVisitationFrequency(relevantSites.at(0));
            particle->setClusterSiteBelongsTo(relevantSites.at(0),clusterId);

          }else{
            mergeSitesToCluster_(relevantSites,favoredClusterId);
            for(int index=1; index < static_cast<int>(relevantSites.size());++index){
              particle->removeMemory(relevantSites.at(index));
            }
            particle->resetVisitationFrequency(relevantSites.at(0));
            particle->setClusterSiteBelongsTo(relevantSites.at(0),favoredClusterId);
          }

        }
      }
    }
  }


}
