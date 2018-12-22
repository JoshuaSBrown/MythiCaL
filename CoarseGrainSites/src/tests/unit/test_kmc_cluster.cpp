#include <iostream>
#include <cassert>
#include <vector>
#include <memory>
#include <cmath>

#include "../../libkmccoarsegrain/topologyfeatures/kmc_cluster.hpp"
#include "../../libkmccoarsegrain/topologyfeatures/kmc_site.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

//  Site site;

  cout << "Testing: Cluster constructor" << endl;
  {
    KMC_Cluster cl;
  }

  cout << "Testing: Cluster identity setter" << endl;
  {
    KMC_Cluster cl;
    cl.setId(0);
  }

  cout << "Testing: Cluster identity getter" << endl;
  {
    KMC_Cluster cl;
    cl.setId(0);
    assert(cl.getId()==0);

    bool fail = false;
    Identity cl2;
    try {
      cl2.getId();
    }catch(...){
      fail = true;
    }
    assert(fail);
  }  

  cout << "Testing: addSite" << endl;
  {
    KMC_Site site;
    site.setId(1);
    double rate = 1.0;
    site.addNeighRate(pair<int, double *>(2,&rate));

    KMC_Cluster cluster;
    cluster.addSite(site);
    assert(cluster.getNumberOfSitesInCluster()==1);
  }


  cout << "Testing: siteIsInCluster" << endl;
  {
    KMC_Site site;
    site.setId(1);
    
    KMC_Cluster cluster;
    assert(!cluster.siteIsInCluster(1));
    cluster.addSite(site);
    assert(cluster.siteIsInCluster(1));
  }

  cout << "Testing: getVisitFrequency" << endl;
  {

    double rate = 1.0;
    KMC_Site site;
    site.setId(1);
    site.setVisitFrequency(5);
    site.addNeighRate(pair<int,double *>(2,&rate));
    assert(site.getVisitFrequency()==5);

    KMC_Site site2;
    site2.setId(2);
    site2.setVisitFrequency(3);
    site2.addNeighRate(pair<int,double *>(1,&rate));
    site2.addNeighRate(pair<int,double *>(3,&rate));
    assert(site2.getVisitFrequency()==3);

    KMC_Cluster cluster;
    assert(!cluster.siteIsInCluster(1));
    cluster.addSite(site);
    cluster.addSite(site2);
    cluster.updateProbabilitiesAndTimeConstant();
    assert(cluster.siteIsInCluster(1));
 
    cout << "Frequency " << cluster.getVisitFrequency(1) << endl;
    assert(cluster.getVisitFrequency(1)==0);
  }
 
  cout << "Testing: setVisitFrequency" << endl;
  {
    double rate = 1.0;
    KMC_Site site;
    site.setId(1);
    site.setVisitFrequency(5);
    site.addNeighRate(pair<int,double *>(2,&rate));
    assert(site.getVisitFrequency()==5);

    KMC_Site site2;
    site2.setId(2);
    site2.setVisitFrequency(3);
    site2.addNeighRate(pair<int,double *>(1,&rate));
    site2.addNeighRate(pair<int,double *>(3,&rate));
    assert(site2.getVisitFrequency()==3);

    KMC_Cluster cluster;
    assert(!cluster.siteIsInCluster(1));
    cluster.addSite(site);
    cluster.addSite(site2);
    cluster.updateProbabilitiesAndTimeConstant();
    assert(cluster.siteIsInCluster(1));
 
    cout << "Frequency " << cluster.getVisitFrequency(1) << endl;
    assert(cluster.getVisitFrequency(1)==0);

    cluster.setVisitFrequency(3,1);
    assert(cluster.getVisitFrequency(1)==3);
    assert(site.getVisitFrequency()==5);
  }   

  cout << "Testing: getProbabilityOfOccupyingInternalSite 1" << endl;
  {

    // Simple convergence test 
    // 
    // site1 -> site2
    //       <-
    //
    // Same rate should lead to 50 % probability on either site

    KMC_Site site;
    site.setId(1);
    double rate = 1;
    site.addNeighRate(pair<int, double *>(2,&rate));
    
    KMC_Site site2;
    site2.setId(2);
    double rate2 = 1;
    site2.addNeighRate(pair<int, double *>(1,&rate2));
  
    KMC_Cluster cluster;
    cluster.addSite(site);
    cluster.addSite(site2);
    cluster.updateProbabilitiesAndTimeConstant();

    assert(static_cast<int>(round(100*cluster.getProbabilityOfOccupyingInternalSite(1)))==50);
    assert(static_cast<int>(round(100*cluster.getProbabilityOfOccupyingInternalSite(2)))==50);

  }

  cout << "Testing: getProbabilityOfOccupyingInternalSite 2" << endl;
  {

    // Simple convergence test 
    // 
    // site1 -> site2  -> site3
    //       <-        <-
    //
    // Same rate should lead to 25 % probability on end sites
    // and 50 % probability on middle site

    KMC_Site site;
    site.setId(1);
    double rate = 1;
    site.addNeighRate(pair<int, double *>(2,&rate));
    
    KMC_Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int, double *>(1,&rate2));
    site2.addNeighRate(pair<int, double *>(3,&rate3));
  
    KMC_Site site3;
    site3.setId(3);
    double rate4 = 1;
    site3.addNeighRate(pair<int, double *>(2,&rate4));

    assert(static_cast<int>(10*site.getProbabilityOfHoppingToNeighboringSite(2))==10);
    assert(static_cast<int>(10*site2.getProbabilityOfHoppingToNeighboringSite(1))==5);
    assert(static_cast<int>(10*site2.getProbabilityOfHoppingToNeighboringSite(3))==5);
    assert(static_cast<int>(10*site3.getProbabilityOfHoppingToNeighboringSite(2))==10);

    KMC_Cluster cluster;
    cluster.setConvergenceIterations(10);

    cluster.addSite(site);
    cluster.addSite(site2);
    cluster.addSite(site3);
    cluster.updateProbabilitiesAndTimeConstant();

    vector<int> values;
    values.push_back(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(1))));
    values.push_back(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(2))));
    values.push_back(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(3))));

    int count_value_25 = 0;
    int count_value_50 = 0;
    for(auto value : values){
      if(value==25) ++count_value_25;
      if(value==50) ++count_value_50;
    }
    cout << "count value 25 " << count_value_25 << endl;
    assert(count_value_25==2);
    assert(count_value_50==1);

  }

  cout << "Testing: get probability of hopping to a neighbor" << endl;
  {

    // Neighbors are not in the cluster
    //
    // neigh1 -> site2  -> site3 -> neigh4
    //        <-        <-       <-
 
    KMC_Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int, double *>(1,&rate2));
    site2.addNeighRate(pair<int, double *>(3,&rate3));
  
    KMC_Site site3;
    site3.setId(3);
    double rate4 = 1;
    double rate5 = 1;
    site3.addNeighRate(pair<int, double *>(2,&rate4));
    site3.addNeighRate(pair<int, double *>(4,&rate5));

    KMC_Cluster cluster;
    cluster.setConvergenceIterations(6);

    cluster.addSite(site2);
    cluster.addSite(site3);

    cluster.updateProbabilitiesAndTimeConstant(); 
    // 
    // site1 -> site2  -> site3 -> neigh4
    //       <-        <-       <-
    assert(static_cast<int>(round(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(1))) ==50);
    assert(static_cast<int>(round(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(4))) ==50);

  }

  cout << "Testing: getProbabilityOfHoppingToNeighborOfCluster 2" << endl;
  {
    KMC_Site site;
    site.setId(1);
    double rate = 1;
    double rate6 = 1;
    site.addNeighRate(pair<int, double *>(2,&rate));
    site.addNeighRate(pair<int, double *>(5,&rate6));
 
    KMC_Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int, double *>(1,&rate2));
    site2.addNeighRate(pair<int , double *>(3,&rate3));
  
    KMC_Site site3;
    site3.setId(3);
    double rate4 = 1;
    site3.addNeighRate(pair<int , double *>(2,&rate4));
    double rate5 = 1;
    site3.addNeighRate(pair<int , double *>(4,&rate5));

    KMC_Cluster cluster;
    cluster.setConvergenceIterations(6);

    cluster.addSite(site);
    cluster.addSite(site2);
    cluster.addSite(site3);

    // Setting the seed will ensure that the results are reproducable
    cluster.updateProbabilitiesAndTimeConstant(); 
    // 
    // neigh5 <- site1 -> site2  -> site3 -> neigh4
    //                 <-        <-
    assert(static_cast<int>(round(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(4))) ==50);
    assert(static_cast<int>(round(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(5))) ==50);

  }

  cout << "Testing: pickNewSiteId" << endl;
  {
    KMC_Site site;
    site.setId(1);
    double rate = 1;
    double rate6 = 1;
    site.addNeighRate(pair<int, double *>(2,&rate));
    site.addNeighRate(pair<int, double *>(5,&rate6));
 
    KMC_Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int, double *>(1,&rate2));
    site2.addNeighRate(pair<int , double *>(3,&rate3));
  
    KMC_Site site3;
    site3.setId(3);
    double rate4 = 1;
    site3.addNeighRate(pair<int , double *>(2,&rate4));
    double rate5 = 1;
    site3.addNeighRate(pair<int , double *>(4,&rate5));

    KMC_Cluster cluster;
    cluster.setConvergenceIterations(6);
    cluster.setResolution(2);
    cluster.addSite(site);
    cluster.addSite(site2);
    cluster.addSite(site3);
    cluster.updateProbabilitiesAndTimeConstant();
    // Setting the seed will ensure that the results are reproducable

    cluster.setRandomSeed(1);
   
    int total = 20000;
    int site1chosen = 0;
    int site2chosen = 0;
    int site3chosen = 0;
    int site4chosen = 0;
    int site5chosen = 0;

    int site1_id = 1;
    int site2_id = 2;
    int site3_id = 3;

    int walker_id = 1;
    cluster.getDwellTime(walker_id);
    for(int count=0; count < total; ++count){
      int chosen_site = cluster.pickNewSiteId(walker_id);
      switch(chosen_site)
      {
        case 1: ++site1chosen;
                cluster.occupy(site1_id);
                cluster.getDwellTime(walker_id);
                break;
        case 2: ++site2chosen;
                cluster.occupy(site2_id);
                cluster.getDwellTime(walker_id);
                break;
        case 3: ++site3chosen;
                cluster.occupy(site3_id);
                cluster.getDwellTime(walker_id);
                break;
        case 4: ++site4chosen;
                cluster.getDwellTime(walker_id);
                break;
        case 5: ++site5chosen;
                cluster.getDwellTime(walker_id);
                break;
      }
    }

    float percent1 = static_cast<float>(site1chosen)/static_cast<float>(total);
    float percent2 = static_cast<float>(site2chosen)/static_cast<float>(total);
    float percent3 = static_cast<float>(site3chosen)/static_cast<float>(total);
    float percent4 = static_cast<float>(site4chosen)/static_cast<float>(total);
    float percent5 = static_cast<float>(site5chosen)/static_cast<float>(total);

    // Able to reproduce the proportion of time spent on each site
    assert(percent1<0.185 && percent1>0.16);
    assert(percent2<0.255 && percent2>0.13);
    assert(percent3<0.185 && percent3>0.165);
    assert(percent4<0.21 && percent4>0.19);
    assert(percent5<0.205 && percent5>0.18);

    vector<double> visit_freq(5,0.0);
    visit_freq.at(0) = static_cast<double>(cluster.getVisitFrequency(1));
    visit_freq.at(1) = static_cast<double>(cluster.getVisitFrequency(2));
    visit_freq.at(2) = static_cast<double>(cluster.getVisitFrequency(3));

    double sum=0.0;
    for(auto visits : visit_freq){
      sum+=visits;
    }
    vector<double> visit_prob;
    for(auto visits : visit_freq){
      visit_prob.push_back(visits/sum);
    }

    assert(visit_prob.at(0) < 0.31);
    assert(visit_prob.at(0) > 0.28);
    assert(visit_prob.at(1) < 0.42);
    assert(visit_prob.at(1) > 0.40);
    assert(visit_prob.at(2) < 0.31);
    assert(visit_prob.at(2) > 0.28);
  }

	return 0;
}
