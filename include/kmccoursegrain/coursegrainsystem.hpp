#ifndef KMCCOURSEGRAIN_COURSEGRAINSYSTEM_H_
#define KMCCOURSEGRAIN_COURSEGRAINSYSTEM_H_

#include <memory>
#include <map>

namespace kmccoursegrain{

class Site;

class CourseGrainSystem{

  public:
    CourseGrainSystem() : tolerance_(0.01) {};
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
        std::map<int,std::map<int,double *>> rates_of_all_sites);

    void visitSite(Particle & particle,int siteId);

    void setCourseGrainingThreshold(int threshold);
    // This will update the internals when the rates vary by the set tolerance 
    void setTheUpdateRateToleranceThreshold(double tolerance);
  private:

    // How many times does a charge have to visit the same sites before it should be course grained
    int courseGrainingThreshold_;
    double tolerance_;
    std::map<int,std::unique_ptr<Site>> sites_;    
    std::map<int,std::unique_ptr<Cluster>> clusters_;
};

}
#endif // KMCCOURSEGRAIN_COURSEGRAINSYSTEM_H_
