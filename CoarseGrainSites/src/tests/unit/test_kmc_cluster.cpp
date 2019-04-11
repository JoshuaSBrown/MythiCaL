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

    double value = static_cast<int>(round(100*cluster.getProbabilityOfOccupyingInternalSite(1)));
    assert(value==50);
    value = static_cast<int>(round(100*cluster.getProbabilityOfOccupyingInternalSite(2)));
    assert(value==50);

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
    double rate6 = 0.001;
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
    double rate5 = 0.001;
    site3.addNeighRate(pair<int , double *>(4,&rate5));

    KMC_Site site4;
    site4.setId(4);
    double rate7 = 1;
    site4.addNeighRate(pair<int, double *>(3,&rate7));

    KMC_Site site5;
    site5.setId(5);
    double rate8 = 1;
    site5.addNeighRate(pair<int, double *>(1,&rate8));

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
    int site4_id = 4;
    int site5_id = 5;

    // 
    // neigh5 -> site1 -> site2  -> site3 -> neigh4
    //        <-       <-        <-       <-
    int walker_id = 1;
    { // Testing cluster pickNewSiteId
      int chosen_site = 2;
      for(int count=0; count < total; ++count){
        switch(chosen_site)
        {
          case 1: ++site1chosen;
                  cluster.occupy(site1_id);
                  cluster.getDwellTime(walker_id);
                  chosen_site = cluster.pickNewSiteId(walker_id);
                  cluster.vacate(walker_id);
                  break;
          case 2: ++site2chosen;
                  cluster.occupy(site2_id);
                  cluster.getDwellTime(walker_id);
                  chosen_site = cluster.pickNewSiteId(walker_id);
                  cluster.vacate(walker_id);
                  break;
          case 3: ++site3chosen;
                  cluster.occupy(site3_id);
                  cluster.getDwellTime(walker_id);
                  chosen_site = cluster.pickNewSiteId(walker_id);
                  cluster.vacate(walker_id);
                  break;
          case 4: ++site4chosen;
                  site4.occupy(site4_id);
                  chosen_site = site4.pickNewSiteId(); 
                  site4.vacate(site4_id);
                  break;
          case 5: ++site5chosen;
                  site5.occupy(site5_id);
                  chosen_site = site5.pickNewSiteId(); 
                  site5.vacate(site5_id);
                  break;
        }
      }
    }

    int baseline_site1chosen = 0;
    int baseline_site2chosen = 0;
    int baseline_site3chosen = 0;
    int baseline_site4chosen = 0;
    int baseline_site5chosen = 0;
    { // Using KMC site pick New site as baseline
      int chosen_site = 2;
      for(int count=0; count < total; ++count){
        switch(chosen_site)
        {
          case 1: ++baseline_site1chosen;
                  site.occupy(site1_id);
                  site.getDwellTime(walker_id);
                  chosen_site = site.pickNewSiteId(); 
                  site.vacate(site1_id);
                  break;
          case 2: ++baseline_site2chosen;
                  site2.occupy(site2_id);
                  site2.getDwellTime(walker_id);
                  chosen_site = site2.pickNewSiteId(); 
                  site2.vacate(site2_id);
                  break;
          case 3: ++baseline_site3chosen;
                  site3.occupy(site3_id);
                  site3.getDwellTime(walker_id);
                  chosen_site = site3.pickNewSiteId(); 
                  site3.vacate(site3_id);
                  break;
          case 4: ++baseline_site4chosen;
                  site4.occupy(site4_id);
                  site4.getDwellTime(walker_id);
                  chosen_site = site4.pickNewSiteId(); 
                  site4.vacate(site4_id);
                  break;
          case 5: ++baseline_site5chosen;
                  site5.occupy(site5_id);
                  site5.getDwellTime(walker_id);
                  chosen_site = site5.pickNewSiteId(); 
                  site5.vacate(site5_id);
                  break;
        }
      }
    }

    float baseline_percent1 = static_cast<float>(site1chosen)/static_cast<float>(total);
    float baseline_percent2 = static_cast<float>(site2chosen)/static_cast<float>(total);
    float baseline_percent3 = static_cast<float>(site3chosen)/static_cast<float>(total);
    float baseline_percent4 = static_cast<float>(site4chosen)/static_cast<float>(total);
    float baseline_percent5 = static_cast<float>(site5chosen)/static_cast<float>(total);

    float percent1 = static_cast<float>(site1chosen)/static_cast<float>(total);
    float percent2 = static_cast<float>(site2chosen)/static_cast<float>(total);
    float percent3 = static_cast<float>(site3chosen)/static_cast<float>(total);
    float percent4 = static_cast<float>(site4chosen)/static_cast<float>(total);
    float percent5 = static_cast<float>(site5chosen)/static_cast<float>(total);

    cout << "baseline perc: " << baseline_percent1 << " perc: " << percent1 << endl;
    cout << "baseline perc: " << baseline_percent2 << " perc: " << percent2 << endl;
    cout << "baseline perc: " << baseline_percent3 << " perc: " << percent3 << endl;
    cout << "baseline perc: " << baseline_percent4 << " perc: " << percent4 << endl;
    cout << "baseline perc: " << baseline_percent5 << " perc: " << percent5 << endl;
    
    
    
    // Able to reproduce the proportion of time spent on each site
    assert(percent1< baseline_percent1*1.2 && percent1>baseline_percent1*0.8 );
    assert(percent2< baseline_percent2*1.2 && percent2>baseline_percent2*0.8 );
    assert(percent3< baseline_percent3*1.2 && percent3>baseline_percent3*0.8 );
    assert(percent4< baseline_percent4*1.2 && percent4>baseline_percent4*0.8 );
    assert(percent5< baseline_percent5*1.2 && percent5>baseline_percent5*0.8 ); 

    vector<double> visit_freq(5,0.0);
    visit_freq.at(0) = static_cast<double>(cluster.getVisitFrequency(1));
    visit_freq.at(1) = static_cast<double>(cluster.getVisitFrequency(2));
    visit_freq.at(2) = static_cast<double>(cluster.getVisitFrequency(3));

    vector<double> baseline_visit_freq(5,0.0);
    baseline_visit_freq.at(0) = static_cast<double>(site.getVisitFrequency());
    baseline_visit_freq.at(1) = static_cast<double>(site2.getVisitFrequency());
    baseline_visit_freq.at(2) = static_cast<double>(site3.getVisitFrequency());

    double sum=0.0;
    double baseline_sum = 0.0;
    for(size_t index =0 ; index < visit_freq.size();++index){
      sum+=visit_freq.at(index);
      baseline_sum +=baseline_visit_freq.at(index);
    }
    vector<double> visit_prob;
    vector<double> baseline_visit_prob;
    for(size_t index = 0; index < visit_freq.size();++index){
      visit_prob.push_back(visit_freq.at(index)/sum);
      baseline_visit_prob.push_back(baseline_visit_freq.at(index)/baseline_sum);
      cout << "Baseline visit prob " << baseline_visit_prob.at(index)  
           << " Visit prob: " << visit_prob.at(index) << endl;
    }

    assert(visit_prob.at(0) < baseline_visit_prob.at(0)*1.2);
    assert(visit_prob.at(0) > baseline_visit_prob.at(0)*0.8);
    assert(visit_prob.at(1) < baseline_visit_prob.at(1)*1.2);
    assert(visit_prob.at(1) > baseline_visit_prob.at(1)*0.8);
    assert(visit_prob.at(2) < baseline_visit_prob.at(2)*1.2);
    assert(visit_prob.at(2) > baseline_visit_prob.at(2)*0.8);
  }

	return 0;
}
