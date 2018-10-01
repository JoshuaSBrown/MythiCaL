#ifndef KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP
#define KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP

#include <memory>
#include <map>

namespace kmccoursegrain{

class Site;
class Cluster;
class Particle;

typedef std::shared_ptr<Site> SitePtr;
typedef std::shared_ptr<Cluster> ClusterPtr;

class CourseGrainSystem{

  public:
    CourseGrainSystem() {};
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
        std::map<int,std::map<int const,double * >> ratesOfAllSites);

    
    void hop(Particle & particle);

    void setCourseGrainThreshold(int threshold);
  private:

    // How many times does a charge have to visit the same sites before it should be course grained
    int courseGrainingThreshold_;
    std::map<int,SitePtr> sites_;    


    std::map<int,ClusterPtr> clusters_;

    void courseGrainSiteIfNeeded_(Particle & particle,const int siteId);
    
    void mergeSiteToCluster_(const int siteId, const int clusterId);
    void createCluster_(const int siteId1,const int siteId2);
    void mergeClusters_(const int clusterId1, const int clusterId2);

};

}
#endif // KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP
