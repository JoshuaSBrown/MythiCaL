#include <stdexcept>
#include <iostream>

#include <kmccoursegrain/coursegrainsystem.hpp>
#include <kmccoursegrain/particle.hpp>
#include <kmccoursegrain/site.hpp>

using namespace std;
using namespace kmccoursegrain;

void CourseGrainSystem::initializeSystem(
    map<int,map<int const,double&>> ratesOfAllSites){

  for(auto it = ratesOfAllSites.begin();it!=ratesOfAllSites.end();++it){
    auto site = unique_ptr<Site>(new Site);
    site->setId(it->first);
    site->setRatesToNeighbors(it->second);
    sites_[it->first] = move(site);
  }
}

void CourseGrainSystem::setCourseGrainThreshold(int threshold){
  courseGrainingThreshold_ = threshold;
}

void CourseGrainSystem::setTheUpdateRateToleranceThreshold(double tolerance){
  if(tolerance<0) throw invalid_argument("Tolerance must be positive");
  tolerance_ = tolerance;
}

void CourseGrainSystem::visitSite(Particle & particle,int siteId){

  particle.occupySite(siteId);
  auto frequency = particle.getVisitationFrequencyOfCurrentlyOccupiedSite();

  // If it has exceeded the threshold
  if(frequency > courseGrainingThreshold_){
    cout << "Not yet implemented visitSite for Course grain system" << endl;
    
  }
}
