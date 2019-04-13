#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "../../../include/kmccoarsegrain/kmc_constants.hpp"
#include "../../../include/kmccoarsegrain/kmc_coarsegrainsystem.hpp"
#include "../../../include/kmccoarsegrain/kmc_walker.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

  cout << "Testing: mergeCluster" << endl;
  {
    // In this example we will define 12 sites
    //
    // site1 - site2 - site3 - site4
    //   |       |       |       |
    // site5 - site6 - site7 - site8
    //   |       |       |       |
    // site9 - site10- site11- site12
    //   |       |       |       |
    // site13- site14- site15- site16
    // 
    // Here rates between 6 and 10 are super fast
    //            between 7 and 11 are super fast
    //          
    //            between 6 and 7 are moderate
    //            between 10 and 11 are moderate
    //
    //            all other rates are slow
    // 
    // Two clusters should be found they should then be merged
    //
    double rate_fast = 10000;
    double rate_moderate = 100;
    double rate_slow = 1;

    vector<pair<int,double>> rates;
    rates.push_back(pair<int,double>(2,rate_slow)); // site1->site2
    rates.push_back(pair<int,double>(5,rate_slow)); // site1->site5

    rates.push_back(pair<int,double>(1,rate_slow)); // site2->site1
    rates.push_back(pair<int,double>(3,rate_slow)); // site2->site3
    rates.push_back(pair<int,double>(6,rate_slow)); // site2->site6

    rates.push_back(pair<int,double>(2,rate_slow)); // site3->site2
    rates.push_back(pair<int,double>(4,rate_slow)); // site3->site4
    rates.push_back(pair<int,double>(7,rate_slow)); // site3->site7

    rates.push_back(pair<int,double>(3,rate_slow)); // site4->site3
    rates.push_back(pair<int,double>(8,rate_slow)); // site4->site8
    
    rates.push_back(pair<int,double>(1,rate_slow)); // site5->site1
    rates.push_back(pair<int,double>(6,rate_slow)); // site5->site6
    rates.push_back(pair<int,double>(9,rate_slow)); // site5->site9

    rates.push_back(pair<int,double>(2,rate_slow)); // site6->site2
    rates.push_back(pair<int,double>(5,rate_slow)); // site6->site5
    rates.push_back(pair<int,double>(7,rate_moderate)); // site6->site7 
    rates.push_back(pair<int,double>(10,rate_fast)); // site6->site10

    rates.push_back(pair<int,double>(3,rate_slow)); // site7->site3
    rates.push_back(pair<int,double>(6,rate_moderate)); // site7->site6
    rates.push_back(pair<int,double>(8,rate_slow)); // site7->site8
    rates.push_back(pair<int,double>(11,rate_fast)); // site7->site11

    rates.push_back(pair<int,double>(4,rate_slow)); // site8->site4
    rates.push_back(pair<int,double>(7,rate_slow)); // site8->site7
    rates.push_back(pair<int,double>(12,rate_slow)); // site8->site12

    rates.push_back(pair<int,double>(5,rate_slow)); // site9->site5
    rates.push_back(pair<int,double>(10,rate_slow)); // site9->site10
    rates.push_back(pair<int,double>(13,rate_slow)); // site9->site13

    rates.push_back(pair<int,double>(6,rate_fast)); // site10->site6
    rates.push_back(pair<int,double>(9,rate_slow)); // site10->site9
    rates.push_back(pair<int,double>(11,rate_moderate)); // site10->site11
    rates.push_back(pair<int,double>(14,rate_slow)); // site10->site14

    rates.push_back(pair<int,double>(7,rate_fast)); // site11->site7
    rates.push_back(pair<int,double>(10,rate_moderate)); // site11->site10
    rates.push_back(pair<int,double>(12,rate_slow)); // site11->site12
    rates.push_back(pair<int,double>(15,rate_slow)); // site11->site15

    rates.push_back(pair<int,double>(8,rate_slow)); // site12->site8
    rates.push_back(pair<int,double>(11,rate_slow)); // site12->site11
    rates.push_back(pair<int,double>(16,rate_slow)); // site12->site16

    rates.push_back(pair<int,double>(9,rate_slow)); // site13->site9
    rates.push_back(pair<int,double>(14,rate_slow)); // site13->site14

    rates.push_back(pair<int,double>(10,rate_slow)); // site14->site10
    rates.push_back(pair<int,double>(13,rate_slow)); // site14->site13
    rates.push_back(pair<int,double>(15,rate_slow)); // site14->site15

    rates.push_back(pair<int,double>(11,rate_slow)); // site15->site11
    rates.push_back(pair<int,double>(14,rate_slow)); // site15->site14
    rates.push_back(pair<int,double>(16,rate_slow)); // site15->site16

    rates.push_back(pair<int,double>(12,rate_slow)); // site16->site12
    rates.push_back(pair<int,double>(15,rate_slow)); // site16->site15

    // Ids of each of the sites
    vector<int> idsOfEachSite = { 
      1,  2,  3,  4,
      5,  6,  7,  8,
      9, 10, 11, 12,
     13, 14, 15, 16};

    // Number of neighbors each site has
    vector<int> numberOfNeighbors = {
      2, 3, 3, 2,
      3, 4, 4, 3,
      3, 4, 4, 3,
      2, 3, 3, 2};

    // Each connection between two sites is composed of a rate going to and
    // from the site e.g.
    //
    // site1 - site2 
    //
    // is equivalent to
    //
    // site1 -> site2
    //
    //        +
    //
    // site1 <- site2
    //
    
    // I have chosen to store the rates in a single vector however, they could
    // be stored in any container or multiple objects. 
    
    // Now we are going to store pointers to the doubles in maps  
    unordered_map< int,unordered_map< int,double>> ratesToNeighbors;

    int global_index = 0;
    int number_of_sites = 16;
    for( int index=0;index<number_of_sites;++index){

      unordered_map< int,double> ratesFromSiteToNeighbors;

      for( int rate_index = 0; 
          rate_index<numberOfNeighbors.at(index);
          ++rate_index){

        ratesFromSiteToNeighbors[rates.at(global_index).first] = \
           (rates.at(global_index).second);

        ++global_index;        
      }
      ratesToNeighbors[idsOfEachSite.at(index)] = ratesFromSiteToNeighbors;
    }

    cout << "Running without cluster" << endl;
    double time_limit = 1000;
    // Without cluster
    {
      KMC_CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(1);
      double time_resolution = time_limit/10.0;
      CGsystem.setTimeResolution(time_resolution);
      CGsystem.setMinCoarseGrainIterationThreshold(constants::inf_iterations);
      CGsystem.initializeSystem(ratesToNeighbors);

      class Electron : public KMC_Walker {};

      Electron electron;
      // Place the electron on site 1
      electron.occupySite(1);

      vector<std::pair<int,KMC_Walker>> electrons;
      electrons.push_back(std::pair<int,KMC_Walker>(1,electron));

      CGsystem.initializeWalkers(electrons);

      vector<double> time_spent_on_sites(number_of_sites,0.0);
      vector<int> hops_made_to_sites(number_of_sites,0);

      double time = 0.0;
      int hop_count = 0;
      auto & electron1 = electrons.at(0).second;
      int id = electrons.at(0).first;
      while(time<time_limit){
        CGsystem.hop(id,electron1);
        time_spent_on_sites.at(electron1.getIdOfSiteCurrentlyOccupying()-1) =
          electron1.getDwellTime();
        hops_made_to_sites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)++;
        time +=electron1.getDwellTime();
        ++hop_count;
      }    

      auto clusters = CGsystem.getClusters();
      // There should be no clusters because the threshold has been set so high
      assert(clusters.size()==0);

      double sum_times = 0.0;
      for(auto time_site : time_spent_on_sites) sum_times+=time_site;

      int siteId=1;

      double sum_time_ratio = 0.0;
      vector<double> time_ratios;
      for(auto time_site : time_spent_on_sites){
        time_ratios.push_back(time_site/sum_times);
        ++siteId;
        sum_time_ratio +=time_site/sum_times;
      }

    } // Without cluster formation


    cout << "Running with Cluster" << endl;
    // With cluster formation
    {
      KMC_CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(1);
      double time_resolution = time_limit/10.0;
      CGsystem.setTimeResolution(time_resolution);
      CGsystem.setMinCoarseGrainIterationThreshold(500);
      CGsystem.initializeSystem(ratesToNeighbors);

      class Electron : public KMC_Walker {};

      Electron electron;
      // Place the electron on site 1
      electron.occupySite(1);

      vector<pair<int,KMC_Walker>> electrons;
      electrons.push_back(pair<int,KMC_Walker>(1,electron));

      CGsystem.initializeWalkers(electrons);

      vector<double> time_spent_on_sites(number_of_sites,0.0);
      vector<int> hops_made_to_sites(number_of_sites,0);

      double time = 0.0;
      int hop_count = 0;
      KMC_Walker& electron1 = electrons.at(0).second;
      int id = electrons.at(0).first;
      while(time<time_limit){
        CGsystem.hop(id,electron1);
        time_spent_on_sites.at(electron1.getIdOfSiteCurrentlyOccupying()-1) =
          electron.getDwellTime();
        hops_made_to_sites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)++;
        time += electron1.getDwellTime();
        ++hop_count;
      }    

      vector<vector<int>> clusters = CGsystem.getClusters();
      assert(clusters.size()==1);
      assert(clusters.at(0).size()==4);
      bool site6_found = false;
      bool site7_found = false;
      bool site10_found = false;
      bool site11_found = false;
      for( auto siteId : clusters.at(0) ){
        if(siteId==6) site6_found = true;
        if(siteId==7) site7_found = true;
        if(siteId==10) site10_found = true;
        if(siteId==11) site11_found = true;
      }
      // Check that the appropriate sites have been found in the cluster
      assert(site6_found);
      assert(site7_found);
      assert(site10_found);
      assert(site11_found);

    } // With cluster formation
  }


	return 0;
}
