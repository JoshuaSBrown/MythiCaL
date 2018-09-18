#ifndef KMCCOURSEGRAIN_CLUSTER_H_
#define KMCCOURSEGRAIN_CLUSTER_H_

//#include <iostream>
#include <vector>
#include <map>
#include <memory>
//#include <math.h>

#include <list>

#include "identity.hpp"


//All functions return 1 if success, 0 or nan if failure as well as cerr printing
namespace kmccoursegrain{

  class Site;

  class Cluster : public virtual Identity {
    public:	
	
      Cluster();

      ~Cluster() {};
      void addSite(std::shared_ptr<Site> siteToAdd);

      bool siteIsInCluster(int siteId) { return sitesInCluster_.count(siteId);}
     
      std::vector<std::shared_ptr<Site>> getSitesInCluster();
      //Overloaded function from Site, prints cluster info
      int getNumberOfSitesInCluster() { return sitesInCluster_.size(); }
      double dwellTime();
//      int convergence(long iterations);
      //Two ways of calling probHopOff, target is target site it leaves the cluster
//      double probHopOff(std::shared_ptr<Site> target, long interations);
//      double probHopOff(int targetId, long interations);
      friend std::ostream& operator<<(std::ostream& os, const kmccoursegrain::Cluster& cluster);
    private:
      std::map<const int, std::shared_ptr<Site>> sitesInCluster_;
//      std::list<int> sitesInCluster_;
      int visitFreqCluster_;


      std::map<const int,double> probabilityOnSite_;

      // first int is hopping from
      // second int is id hopping too
      // double is the rate
      std::map<const int,std::vector<std::pair<const int, double>>> getInternalRates();

      // First int is the Id of a site within the cluster
      // pair - first int is the id of the site neighboring the cluster
      // double is the rate
//      std::map<int,std::pair<int,double>> getRatesToNeighborsOfCluster_();
 
      void initializeProbabilityOnSite_(); 

      //Notes:
      //Stuff to do not in matlab
      //smoosh sites together if below thresh
      //store probilties and data somewhere
      //Prob hops to a ceritain neighbor site off cluster, using probality from cluster convergence

      //calculate site ratio given hop rates off site, need list of neighbors for that site, hop rates to the neigh form site, list of sites Ids in cluster
      //pass maybe cluster struct and return array of site ratios of the sites hop off/hop on see matlab file
      //
      //fn 2 determine sites in cluster and return list of id sites in cluster,
      //
      //struct cluster of 2 lsits ids in cluster and neighbors, int id of cluster, int # of sites in cluster, int # of neighbors
      //
      //create a repository and error handeling
      //
      //dwell time 1/sum of rates; list of neighbor rates; return double
      //
      //calculate pvals for site rate(#)/sum of rates site #; list of neighbor hop rates, out put an array
      //do for all sites in and out the cluster
      //
      //calculate p hop off given lisst of neigh rates for a site, output an array
      //
      //
      //Look here 
      //dwell time = 1/(sum of rates to neighbors)
      //site prob hop = rate to specific site / (sum of rates to all neighbors)
      //prob to hop off cluster = rate to hop to neighbor from cluster / (sum of rates to all neighbors to cluster)
      //NOTE: prob hop off cluster to site in cluster is 0
      //prob on site INITIALLY is 1/number of sites in cluster
      //run convergence to fix prob on site in cluster
      //t escape = 1/(rates to escape from site out of cluster)
      //prob off cluster on neighbor = prob on site in cluster * (dwell time of site in cluster/total dwell time of all sites in cluster) * (rate off to neighbor from site/ total rate off cluster to neighbors)
      //
      //Return prob on site in cluster prob off cluster on neighbor
      //
      //convergence
      //	1. prob on site is 1/(number of sites in cluster) (only intially)
      //	2. site 1 eg (prob hope to site 1)*(prob on site 2)+(prob hppe to site 1)(prob on site 5)
      //	3.Total (new site prob)     total=site1prob+site2prob...
      //	3. normalize
      //	4. Average (site1oldProb + site1newProb)/2
      //	4.probsite1=site1prob/total + probsite1 ....
      //	5.normalise probonsite1=probonsite1/sum(probs on site #)
      //
      //calculate total prob off cluster to a given neighbor n
      //	1. sum hop rates off the cluster
      //	2. calc prob hop to neigh off cluster ex for neigh1 prob neigh=dwellofsite(1)/sum(dwell + probonsite1*(rate from site1 to neigh1/tot    tot =sum of hop rates off
      //
      //is site(1) above refer to as the cluster or a sight in the cluster?
      //
      //note: if charge is within cluster, most likeyl on which cluster, where most likely to jump off cluster to which site?, dwell time on cluster as single site
      //
      //
      //prob to hop off from a given site in cluster n
      //hopoff1 = ProbOnSite1 * Dwell1/sum(dwell) * rate from 1to neighbor n + probonsite1*dwell/sum(dwell)*rate from site1toneigh m
      //then normalize hopoff[i]/sum(hopoff)
      //
      //escape time off cluster
      //	1. prob hop off site i time =1/(rate to neigh n + rate to neigh m)
      //	escape time = site1hop off time * hop off site 1+...
      //
      //Tprob_off cluster or escape cluster = tescapsesite1*probhopoffsite1+...
  };

  // sets thresh as a static for all functions, not in class cluster
  void setThreshold(const int n);

  //returns the threshold
  int getThreshold();

//  static bool siteAboveThreshold(int siteId, int frequency_visitation);

//  Cluster generateCluster(int clusterId);

}

#endif // KMCCOURSEGRAIN_CLUSTER_H
