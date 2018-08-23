#include <iostream>
#include <kmccoursegrain/coursegrainsystem.hpp>

using namespace std;

using namespace kmccoursegrain;

void CourseGrainSystem::initializeSystem(
    map<int,map<int,double*>> rates_of_all_sites){

  for(auto map_pr : rates_of_all_sites){
    Site site;
    site.setId(map_pr.first);
    site.setRatesToNeighbors(map_pr.second);
    sites_[map_pr.first] = make_unique(site);
  }
}

void CourseGrainSystem::setCourseGrainThreshold(int threshold){
  courseGrainingThreshold_(threshold);
}

void CourseGrainSystem::setTheUpdateRateToleranceThreshold(double tolerance){
  tolerance_ = tolerance;
}

void CourseGrainSystem::visitSite(Particle & particle,int siteId){

  particle.recordSiteVisitEvent(siteId);
  auto particleStats = particle.getSiteIdWithMaxVisitCount();

  // If it has exceeded the threshold
  if(particleStats.second > courseGrainingThreshold_){
    cout << "Not yet implemented visitSite for Course grain system" << endl;
  }
}
