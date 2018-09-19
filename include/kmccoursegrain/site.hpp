#ifndef KMCCOURSEGRAIN_SITE_H_
#define KMCCOURSEGRAIN_SITE_H_

#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <math.h>

#include "identity.hpp"

namespace kmccoursegrain {
  //All functions return 1 if success, 0 or nan if failure as well as cerr printing
  class Site : public virtual Identity {
    //	friend class cluster;
    public:	
      Site() : visitFreq_(0), clusterId_(-1) {};
      void setRatesToNeighbors(std::map<int const, double*> neighRates);

      void addNeighRate(std::pair<const int, double*> neighRate);

      void resetNeighRate(std::pair<const int,double*> neighRate);

      bool isNeighbor(int neighId) { return neighRates_.count(neighId);}

      std::vector<double> getRatesToNeighbors();
      double getRateToNeighbor(int neighId) { return *(neighRates_[neighId]);}

      std::vector<int> getNeighborIds();

      double getDwellTime();
      //std::vector<int> getIdsOfNeighSites();
      /*
      //Two ways to create a site object
      //Creates site with list of sites and rates as neighbors
      site(int sId, int vFreq, std::map<std::shared_ptr<site>,double> nSites);
      //Just creates site
      site(int sId, int vFreq);*/
      //Two ways to call probHop hop is from shipping to receivingSite
//      double probHop(std::shared_ptr<site> receivingSite);
      void setClusterId(const int clusterId) { clusterId_ = clusterId; }
      double probHopToNeigh(const int neighSiteId);
/*      int addNeighbors(std::map<std::shared_ptr<site>, double>);
      void printInfo();*/
//      int mergeSites();//FIXME

      friend std::ostream& operator<<(std::ostream& os, const kmccoursegrain::Site& site);
    private:
      //		int siteId;
      //		std::vector<std::shared_ptr<site>> neighSites;
      //		std::map<int, double> neighs;
      std::map<int const, double *> neighRates_;
      int visitFreq_;
   //   double courseGrainedSiteProb_;
   //   bool probabilitySet_;
      int clusterId_;
      //		double probOnSite;
  };

}


#endif
