#ifndef KMCCOURSEGRAIN_SITE_H_
#define KMCCOURSEGRAIN_SITE_H_

#include <map>
#include "identity.hpp"

namespace kmccoursegrain{

//class Cluster;

class Site : public Identity {
  public:
    Site() : probabilitySet_(false), visitFreq_(0) {}; //, cluster_ptr_(nullptr) {};
    void setRatesToNeighbors(std::map<int const,double&> neighRates);
//    void setCluster(Cluster * cluster);
  private:

    int visitFreq_;
//    Cluster * cluster_ptr_;
    double courseGrainedSiteProb_;
    bool probabilitySet_;
    std::map<int const,double&> neighRates_;

};

}
#endif // KMCCOURSEGRAIN_SITE_H_
