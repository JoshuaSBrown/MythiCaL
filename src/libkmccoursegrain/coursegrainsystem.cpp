#include <stdexcept>
#include <iostream>

#include "../../include/kmccoursegrain/coursegrainsystem.hpp"
#include "../../include/kmccoursegrain/particle.hpp"
#include "kmc_site.hpp"
#include "log.hpp"
#include "kmc_cluster.hpp"

using namespace std;

namespace kmccoursegrain {

  void CourseGrainSystem::initializeSystem(
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

  void CourseGrainSystem::initializeParticles(vector<Particle*> particles){

    LOG("Initializeing particles",1);

    if(sites_.size()==0){
      throw runtime_error("You must first initialize the system before you "
          "can initialize the particles");
    }

    for(auto particle : particles){
      auto siteId = particle->getIdOfSiteCurrentlyOccupying();
      if(siteId == kmc_particle::unassignedSiteId){
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

  void CourseGrainSystem::setCourseGrainThreshold(int threshold){
    LOG("Setting threshold",1);
    courseGrainingThreshold_ = threshold;
  }

  void CourseGrainSystem::setRandomSeed(const unsigned long seed){
    if(sites_.size()!=0){
      throw runtime_error("For the random seed to have an affect, it must be "
          "set before initializeSystem is called");
    }
    seed_=seed;
    seed_set_=true;
  }

  void CourseGrainSystem::createCluster_(const int siteId1,const int siteId2){
    LOG("Creating cluster",1);
    auto cluster_ptr = shared_ptr<KMC_Cluster>(new KMC_Cluster());
    cluster_ptr->addSite(sites_[siteId1]);
    cluster_ptr->addSite(sites_[siteId2]);
    if(seed_set_){
      cluster_ptr->setRandomSeed(seed_);
      ++seed_;
    }
    clusters_[cluster_ptr->getId()] = cluster_ptr;
  }

  void CourseGrainSystem::mergeSiteToCluster_(
      const int siteId,
      const int clusterId){
    LOG("Merging site to cluster",1);

    clusters_[clusterId]->addSite(sites_[siteId]);
  }

  void CourseGrainSystem::mergeClusters_(
      const int clusterId1,
      const int clusterId2){
    LOG("Merging clusters",1);

    if(clusterId1==clusterId2){
      throw invalid_argument("You cannot merge two clusters with the same Id");
    }

    if(clusterId1<clusterId2){
      clusters_[clusterId1]->migrateSitesFrom(clusters_[clusterId2]);
      clusters_.erase(clusterId2); 
    }else{
      clusters_[clusterId2]->migrateSitesFrom(clusters_[clusterId1]);
      clusters_.erase(clusterId1); 
    }
  }

  void CourseGrainSystem::courseGrainSiteIfNeeded_(Particle * particle){

    LOG("Course graining sites if needed",1);

    auto siteFrequencies = particle->getVisitationFrequenciesOfCurrentSites();

    // Determine if the particle has made enough jumps to have a memory
    if(siteFrequencies.size()>1){
      // Get frequencies of visitiation of the currently occupied site and 
      // the site it just hopped from  
      int siteId1 = siteFrequencies.at(0).first;
      int siteId2 = siteFrequencies.at(1).first;
      int frequency1 = siteFrequencies.at(0).second;
      int frequency2 = siteFrequencies.at(1).second;
      // If it has exceeded the threshold

      if(frequency1 > courseGrainingThreshold_ &&
          frequency2 > courseGrainingThreshold_){

        // Determine if the two sites are already part of a cluster
        int clusterId1 = sites_[siteId1]->getClusterId();
        int clusterId2 = sites_[siteId2]->getClusterId();

        // If neither is part of a cluster create one
        if(clusterId1==-1 && clusterId2==-1){
          createCluster_(siteId1,siteId2);
        }else if(clusterId1==-1 && clusterId2>-1){
          // site 1 is not part of a cluster but site 2 is 
          // join site 1 to cluster of site 2
          mergeSiteToCluster_(siteId1,clusterId2);

        }else if(clusterId1>-1 && clusterId2==-1){
          // site 1 is part of a cluster but site 2 is not
          // join site 2 to cluster of site 1
          mergeSiteToCluster_(siteId2,clusterId1);

        }else if(clusterId1!=clusterId2){
          // both sites are exist on two different clusters
          // merge the cluster with the larger id with the cluster
          // that has the smaller id
          mergeClusters_(clusterId1,clusterId2);
        }

      }
    }
  }

  void CourseGrainSystem::hop(Particle * particle){

    LOG("Particle is hopping in system",1);
    auto siteId = particle->getIdOfSiteCurrentlyOccupying();
    int siteToHopTo = particle->getPotentialSite();

    // Determine if the potential site is occupied or not
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
        courseGrainSiteIfNeeded_(particle);

        int newId = clusters_[clusterId]->pickNewSiteId();
        auto hopTime = clusters_[clusterId]->getDwellTime();

        particle->occupySite(siteToHopTo);
        particle->setDwellTime(hopTime);
        particle->setPotentialSite(newId);

      }else{

        sites_[siteId]->vacateSite();
        sites_[siteToHopTo]->occupySite();
        courseGrainSiteIfNeeded_(particle);

        int newId = sites_[siteToHopTo]->pickNewSiteId();
        auto hopTime = sites_[siteToHopTo]->getDwellTime();

        particle->occupySite(siteToHopTo);
        particle->setDwellTime(hopTime);
        particle->setPotentialSite(newId);
      }

    }

  }

}
