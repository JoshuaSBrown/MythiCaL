#ifndef KMCCOURSEGRAIN_PARTICLE_H_
#define KMCCOURSEGRAIN_PARTICLE_H_

#include <list>

namespace kmccoursegrain{

/**
 * \brief abstract class meant to be inherited 
 */
class Particle{
  public: 
    Particle() {};
    virtual void setMemoryCapacity(size_t numberSitesToRemember);
    virtual size_t getMemoryCapacity();
    virtual void occupySite(int siteId);
    virtual int getVisitationFrequencyOfCurrentlyOccupiedSite();
    virtual std::list<std::pair<int,int>> getMemory();
  private:
    std::list<std::pair<int,int>> memoryQueue_;

    bool rememberSite_(int siteId);
    std::list<std::pair<int,int>>::iterator getMemoryIterator_(int siteId);
    void refreshMemory_(int siteId);
    void createNewMemory_(int siteId);
};

}
#endif // KMCCOURSEGRAIN_PARTICLE_H_
