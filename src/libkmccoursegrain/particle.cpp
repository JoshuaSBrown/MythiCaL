#include <stdexcept>
#include <cassert>

#include "../../include/kmccoursegrain/particle.hpp"

using namespace std;

namespace kmccoursegrain {

  /****************************************************************************
   * External Methods
   ****************************************************************************/

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

  unsigned int Particle::getMemoryCapacity() const {
    return static_cast<unsigned int>(memoryQueue_.size());
  }

  int Particle::getVisitationFrequencyOfCurrentlyOccupiedSite(){
    assert(memoryQueue_.size()>0);
    for( auto memory : memoryQueue_ ){
      if(memory.first != unassignedSiteId_) return memory.second;
    }
    string err = "ERROR you cannot get the visitation frequency of the "
      "unoccupied site at this moment, as no sites have yet been occupied.\n";
    throw runtime_error(err);
  }

  void Particle::occupySite(const int siteId){
    assert(siteId!=unassignedSiteId_);
    if(rememberSite_(siteId)){
      refreshMemory_(siteId);
    }else{
      createNewMemory_(siteId);
    }
  }

  vector<pair<int,int>> 
  Particle::getVisitationFrequenciesOfCurrentSites() const {

    vector<pair<int,int>> frequencies;

    for(auto memory_it = memoryQueue_.begin();
        memory_it != memoryQueue_.end();
        ++memory_it){

      if(memory_it->first != unassignedSiteId_ ){
        frequencies.push_back(pair<int,int>(memory_it->first,memory_it->second));
        ++memory_it;
        if(memory_it->first != unassignedSiteId_ ){
          frequencies.push_back(pair<int,int>(memory_it->first,memory_it->second));
          return frequencies;
        }else{
          string err = "ERROR only a single site was found to have been occupied"
            " a particle needs to have occupied at least two sites before you "
            "can grab the frequencies of the two most recently occupied sites.";
            throw runtime_error(err);
        }
      }
    }
    string err = "ERROR you cannot get the visitation frequency of the "
      "unoccupied site at this moment, as no sites have yet been occupied.\n";
    throw runtime_error(err);
  }

  int Particle::getIdOfSiteCurrentlyOccupying() const {
    return memoryQueue_.begin()->first;
  }

  list<pair<int,int>> Particle::getMemory(){
    return memoryQueue_;
  }

  /****************************************************************************
   * Internal Methods
   ****************************************************************************/

  bool Particle::rememberSite_(int siteId){
    for(auto siteVisit : memoryQueue_){
      if( siteVisit.first == siteId) return true;
    }
    return false;
  }

  void Particle::refreshMemory_(const int siteId){
    auto it = getMemoryIterator_(siteId);
    ++(it->second);
    memoryQueue_.splice(memoryQueue_.begin(),memoryQueue_,it);
  }

  void Particle::createNewMemory_(const int siteId){
    auto it = memoryQueue_.end();
    --it;
    it->first = siteId;
    it->second = 1;
    memoryQueue_.splice(memoryQueue_.begin(), memoryQueue_,it);
  }

  list<pair<int,int>>::iterator Particle::getMemoryIterator_(const int siteId){
    auto it = memoryQueue_.begin();
    while(it!=memoryQueue_.end()){
      if( it->first == siteId)  return it;
      ++it;
    }
    throw runtime_error("Site does not exist in queue/memory");
  }


}
