#ifndef KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP
#define KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP

#include <vector>
#include <memory>
#include <map>

namespace kmccoursegrain{

class Site;
class Cluster;
class Particle;

typedef std::shared_ptr<Site> SitePtr;
typedef std::shared_ptr<Cluster> ClusterPtr;

/**
 * \brief Course Grain System allows abstraction of renormalization of sites
 *
 * This class simply uses pointers to the relevant hop rates. It will simulate a
 * particle hopping through the system of sites. If a large number of compute
 * cycles are expended moving a particle between two low energy sites they will
 * be course grained, or renormalized so that the probabilities and time spent
 * on each site will be the same but the number of compute cycles will be 
 * significantly reduced. 
 **/
class CourseGrainSystem{

  public:
  
    /**
     * \brief Constuctor for course grained system
     *
     * The constructor by default seeds the random number generator based on the
     * current time. Furthermore, course graining of two sites is by default set
     * to a threshold of 20. Meaning a particle must remember moving back and
     * forth between two sites at least 20 times before the sites are course 
     * grained. 
     **/
    CourseGrainSystem() : seed_set_(false), courseGrainingThreshold_(20) {};
    // There should be no reason to update the rates because we are using 
    // pointers to the rates, though this means the doubles must be stored on
    // the heap somewhere. 
    // The map
    // 
    // int1 - id of site i
    // int2 - id of site j
    // double * - pointer to rate going from i->j
    //
    void initializeSystem(
        std::map<const int,std::map<int const,double * >> ratesOfAllSites);

    void initializeParticles(std::vector<Particle *> particles);

    void setRandomSeed(const unsigned long seed);
    void hop(Particle * particle);

    void setCourseGrainThreshold(int threshold);
  private:

    // How many times does a charge have to visit the same sites before it should be course grained
    bool seed_set_;
    unsigned long seed_;

    int courseGrainingThreshold_;
    std::map<int,SitePtr> sites_;    


    std::map<int,ClusterPtr> clusters_;

    void courseGrainSiteIfNeeded_(Particle * particle);
    
    void mergeSiteToCluster_(const int siteId, const int clusterId);
    void createCluster_(const int siteId1,const int siteId2);
    void mergeClusters_(const int clusterId1, const int clusterId2);

};

}
#endif // KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP
