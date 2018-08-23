#ifndef KMCCOURSEGRAIN_PARTICLE_H_
#define KMCCOURSEGRAIN_PARTICLE_H_

#include <limits>
#include <list>

namespace kmccoursegrain{

/**
 * \brief abstract class meant to be inherited 
 */
class Particle{
  public: 
    Particle() {};
    void setMemoryCapacity(unsigned int numberSitesToRemember);
    unsigned int getMemoryCapacity();// {
//      return static_cast<unsigned int>(memoryQueue_.size()); 
//    }
    void occupySite(int siteId);
    int getVisitationFrequencyOfCurrentlyOccupiedSite();
    std::list<std::pair<int,int>> getMemory();
  private:
    std::list<std::pair<int,int>> memoryQueue_;
    const int unassignedSiteId_ = std::numeric_limits<int>::min();

    bool rememberSite_(int siteId);
    std::list<std::pair<int,int>>::iterator getMemoryIterator_(int siteId);
    void refreshMemory_(int siteId);
    void createNewMemory_(int siteId);
};

}
#endif // KMCCOURSEGRAIN_PARTICLE_H_
