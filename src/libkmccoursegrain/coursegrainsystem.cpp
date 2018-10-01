#include <stdexcept>
#include <iostream>

#include "../../include/kmccoursegrain/coursegrainsystem.hpp"
#include "../../include/kmccoursegrain/particle.hpp"
#include "site.hpp"
#include "cluster.hpp"

using namespace std;

namespace kmccoursegrain {

  void CourseGrainSystem::initializeSystem(
    map<int,map<int const,double* >> ratesOfAllSites){

    for(auto it = ratesOfAllSites.begin();it!=ratesOfAllSites.end();++it){
      auto site = shared_ptr<Site>(new Site);
      site->setId(it->first);
      site->setRatesToNeighbors(it->second);
      sites_[it->first] = move(site);
    }
  }

  void CourseGrainSystem::setCourseGrainThreshold(int threshold){
    courseGrainingThreshold_ = threshold;
  }


  void CourseGrainSystem::createCluster_(const int siteId1,const int siteId2){
    auto cluster_ptr = shared_ptr<Cluster>(new Cluster());
    cluster_ptr->addSite(sites_[siteId1]);
    cluster_ptr->addSite(sites_[siteId2]);
    clusters_[cluster_ptr->getId()] = cluster_ptr;
  }

  void CourseGrainSystem::mergeSiteToCluster_(
      const int siteId,
      const int clusterId){

    clusters_[clusterId]->addSite(sites_[siteId]);
  }

  void CourseGrainSystem::mergeClusters_(
      const int clusterId1,
      const int clusterId2){

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

  void CourseGrainSystem::courseGrainSiteIfNeeded_(Particle & particle){

    auto siteFrequencies = particle.getVisitationFrequenciesOfCurrentSites();
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

  void CourseGrainSystem::hop(Particle & particle){

    auto siteId = particle.getIdOfSiteCurrentlyOccupying();
    int siteToHopTo = particle.getPotentialSite();

    // Determine if the potential site is occupied or not
    if( sites_[siteToHopTo]->siteIsOccupied() ){

      if( sites_[siteId]->partOfCluster() ){

        auto clusterId = sites_[siteId]->getClusterId();
        int newId = clusters_[clusterId]->pickNewSiteId();
        auto hopTime = clusters_[clusterId]->getDwellTime();

        particle.setDwellTime(hopTime);
        particle.setPotentialSite(newId);

      }else{

        int newId = sites_[siteId]->pickNewSiteId();
        auto hopTime = sites_[siteId]->getDwellTime();

        particle.setDwellTime(hopTime);
        particle.setPotentialSite(newId);
      }

    }else{

      if( sites_[siteToHopTo]->partOfCluster() ){

        auto clusterId = sites_[siteToHopTo]->getClusterId();
        sites_[siteId]->vacateSite();
        sites_[siteToHopTo]->occupySite();
        courseGrainSiteIfNeeded_(particle);

        int newId = clusters_[clusterId]->pickNewSiteId();
        auto hopTime = clusters_[clusterId]->getDwellTime();

        particle.occupySite(siteToHopTo);
        particle.setDwellTime(hopTime);
        particle.setPotentialSite(newId);

      }else{

        sites_[siteId]->vacateSite();
        sites_[siteToHopTo]->occupySite();
        courseGrainSiteIfNeeded_(particle);

        int newId = sites_[siteToHopTo]->pickNewSiteId();
        auto hopTime = sites_[siteToHopTo]->getDwellTime();

        particle.occupySite(siteToHopTo);
        particle.setDwellTime(hopTime);
        particle.setPotentialSite(newId);
      }

    }

  }

}
