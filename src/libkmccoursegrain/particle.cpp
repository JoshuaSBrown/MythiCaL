#include <stdexcept>
#include <cassert>

#include <kmccoursegrain/particle.hpp>

using namespace std;

using namespace kmccoursegrain;

/******************************************************************************
 * Internal Methods
 ******************************************************************************/

bool Particle::rememberSite_(int siteId){
  for(auto siteVisit : memoryQueue_){
    if( siteVisit.first == siteId) return true;
  }
  return false;
}

list<pair<int,int>>::iterator Particle::getMemoryIterator_(int siteId){
  auto it = memoryQueue_.begin();
  while(it!=memoryQueue_.end()){
    if( it->first == siteId)  return it;
    ++it;
  }
  throw runtime_error("Site does not exist in queue/memory");
}

void Particle::refreshMemory_(int siteId){
  auto it = getMemoryIterator_(siteId);
  ++(it->second);
  memoryQueue_.splice(memoryQueue_.begin(),memoryQueue_,it);
}

void Particle::createNewMemory_(int siteId){
  auto it = memoryQueue_.end();

  it->first = siteId;
  it->second = 1;

  memoryQueue_.splice(memoryQueue_.begin(), memoryQueue_,it);
}

/******************************************************************************
 * External Methods
 ******************************************************************************/

void Particle::setMemoryCapacity(unsigned int numberSitesToRemember){
  
  if(numberSitesToRemember==1){
    string err = "It does not make sense to have a memory of 1 because a "
      "particle cannot hop to the same site it currently occupies";
    throw invalid_argument(err);
  }

  int delta = static_cast<int>(numberSitesToRemember) - 
              static_cast<int>(memoryQueue_.size());

  if(delta>0){
    for(int i=0;i<delta;++i) {
      memoryQueue_.push_back(pair<int,int>(unassignedSiteId_,0));
    }
  }else{
    delta = delta*-1;
    for(int i=0;i<delta;++i) memoryQueue_.pop_back();
  }
}

unsigned int Particle::getMemoryCapacity(){
  return static_cast<unsigned int>(memoryQueue_.size());
}

void Particle::occupySite(int siteId){
  assert(siteId!=unassignedSiteId_);
  if(rememberSite_(siteId)){
    refreshMemory_(siteId);
  }else{
    createNewMemory_(siteId);
  }
}

int Particle::getVisitationFrequencyOfCurrentlyOccupiedSite(){
  assert(memoryQueue_.size()>0);
  return memoryQueue_.begin()->second;
}

list<pair<int,int>> Particle::getMemory(){
  return memoryQueue_;
}
