#include <stdexcept>
#include <cassert>
#include "../../include/kmccoursegrain/particle.hpp"

using namespace std;

namespace kmccoursegrain {

  /****************************************************************************
   * External Methods
   ****************************************************************************/
  Particle::Particle(){
    potentialSite_ = constants::unassignedId;
    dwelltime_ = -1.0;
    setMemoryCapacity(2);
  }

  void Particle::setMemoryCapacity(unsigned int numberSitesToRemember){

    if(numberSitesToRemember==static_cast<unsigned int>(memoryQueue_.size())) {
      return;
    }

    if(numberSitesToRemember==1){
      string err = "It does not make sense to have a memory of 1 because a "
        "particle cannot hop to the same site it currently occupies";
      throw invalid_argument(err);
    }

    int delta = static_cast<int>(numberSitesToRemember) - 
      static_cast<int>(memoryQueue_.size());

    if(delta>0){
      for(int i=0;i<delta;++i) {
        Memory memory = {
          constants::unassignedId,
          constants::unassignedId,
          0};

        memoryQueue_.push_back(memory);
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
      if(memory.siteId != constants::unassignedId){
        return memory.visitFrequency;
      }
    }
    string err = "ERROR you cannot get the visitation frequency of the "
      "unoccupied site at this moment, as no sites have yet been occupied.\n";
    throw runtime_error(err);
  }

  void Particle::occupySite(
      const int siteId,
      const int clusterId){

    assert(siteId!=constants::unassignedId);

    if(remember_(siteId,clusterId)){
      // Remember both
      refreshMemoryOfCluster_(clusterId);
    }else if(rememberSite_(siteId)){
      // Remember site
      refreshMemoryOfSite_(siteId);
    }else if(rememberCluster_(clusterId)){
      // Remember cluster
      refreshMemoryOfCluster_(clusterId);
      memoryQueue_.begin()->siteId = siteId;
    }else{
      // Remember neither
      createNewMemory_(siteId,clusterId);
    }

  }

  void Particle::setClusterSiteBelongsTo(const int siteId,const int clusterId){
    auto memory_it = getMemoryIteratorSite_(siteId);
    memory_it->clusterId = clusterId;
  }

  vector<vector<int>> Particle::getVisitationFrequenciesOfCurrentSites() const {

    vector<vector<int>> frequencies;
    for(auto memory : memoryQueue_){
      if( memory.siteId != constants::unassignedId ){
        frequencies.push_back(vector<int>{
            memory.siteId, 
            memory.clusterId,
            memory.visitFrequency });
      }
    }

    return frequencies;
  }

  void Particle::resetVisitationFrequency(const int siteId){
    for(auto memory : memoryQueue_){
      if(memory.siteId == siteId ) {
        memory.visitFrequency = 0;
        return;
      }
    }
    throw invalid_argument("Site id is not in the particles memory");
  }

  void Particle::removeMemory(const int siteId){
    auto memory_it = getMemoryIteratorSite_(siteId);
    memory_it->siteId = constants::unassignedId;
    memory_it->clusterId = constants::unassignedId;
    memory_it->visitFrequency = 0;
    memoryQueue_.splice(memoryQueue_.end(),memoryQueue_,memory_it);
  }

  int Particle::getIdOfSiteCurrentlyOccupying() const {
    return memoryQueue_.begin()->siteId;
  }

  list<vector<int>> Particle::getMemory(){
    list<vector<int>> memories;
    for(auto memory : memoryQueue_){
      memories.push_back(vector<int>{
          memory.siteId,
          memory.clusterId,
          memory.visitFrequency});
    }
    return memories;
  }

  int Particle::getPotentialSite() const {
    if(potentialSite_==constants::unassignedId){
      throw runtime_error("Cannot get potential site as it has not yet been "
          "assigned. You many need to first initialize the particle");
    }
    return potentialSite_;
  }
  /****************************************************************************
   * Internal Methods
   ****************************************************************************/

  bool Particle::rememberSite_(const int siteId){

    if(siteId!=constants::unassignedId){
      for(auto memory : memoryQueue_){
        if( memory.siteId == siteId) return true;
      }
    }
    
    return false;
  }

  bool Particle::rememberCluster_(const int clusterId){

    if(clusterId!=constants::unassignedId){
      for(auto memory : memoryQueue_){
        if( memory.clusterId == clusterId) return true;
      }
    }
    return false;
  }

  bool Particle::remember_(const int siteId,const int clusterId){
    if( rememberSite_(siteId) && rememberCluster_(clusterId) ) return true;
    return false;
  }

  void Particle::refreshMemoryOfSite_(const int siteId){
    auto memory_it = getMemoryIteratorSite_(siteId);
    ++(memory_it->visitFrequency);
    memoryQueue_.splice(memoryQueue_.begin(),memoryQueue_,memory_it);
  }

  void Particle::refreshMemoryOfCluster_(const int clusterId){
    auto memory_it = getMemoryIteratorCluster_(clusterId);
    ++(memory_it->visitFrequency);
    memoryQueue_.splice(memoryQueue_.begin(),memoryQueue_,memory_it);
  }

  void Particle::createNewMemory_(const int siteId,const int clusterId){
    auto memory_it = memoryQueue_.end();
    --memory_it;
    memory_it->siteId = siteId;
    memory_it->clusterId = clusterId;
    memory_it->visitFrequency = 1;
    memoryQueue_.splice(memoryQueue_.begin(), memoryQueue_,memory_it);
  }

  list<Particle::Memory>::iterator 
  Particle::getMemoryIteratorSite_(const int siteId){

    auto memory_it = memoryQueue_.begin();
    while(memory_it!=memoryQueue_.end()){
      if( memory_it->siteId == siteId)  return memory_it;
      ++memory_it;
    }
    throw runtime_error("Site does not exist in queue/memory");
  }

  list<Particle::Memory>::iterator 
  Particle::getMemoryIteratorCluster_(const int clusterId){
    auto memory_it = memoryQueue_.begin();
    while(memory_it!=memoryQueue_.end()){
      if( memory_it->clusterId == clusterId)  return memory_it;
      ++memory_it;
    }
    throw runtime_error("Cluster does not exist in queue/memory");
  }


}
