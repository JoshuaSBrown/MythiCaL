#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "mythical/constants.hpp"
#include "mythical/coarsegrainsystem.hpp"
#include "mythical/walker.hpp"

using namespace std;
using namespace mythical;

int main(void){

  cout << "Testing: CoarseGrainSystem constructor" << endl;
  {
    CoarseGrainSystem CGsystem;
  }

  cout << "Testing: initializeSystem" << endl;
  {

    // In this example we will define 12 sites
    //
    // site1 - site2 - site3 - site4
    //   |       |       |       |
    // site5 - site6 - site7 - site8
    //   |       |       |       |
    // site9 - site10- site11- site12
    
    // Define a vector with all the rates starting from site 1 and describing
    // the rates to its neighbors and then continue with the other sites
    
    double rate_fast = 100;
    double rate_slow = 1;

    vector<pair<int,double>> rates;
    rates.push_back(pair<int,double>(2,rate_slow)); // site1->site2
    rates.push_back(pair<int,double>(5,rate_slow)); // site1->site5

    rates.push_back(pair<int,double>(1,rate_slow)); // site2->site1
    rates.push_back(pair<int,double>(3,rate_slow)); // site2->site3
    rates.push_back(pair<int,double>(6,rate_fast)); // site2->site6

    rates.push_back(pair<int,double>(2,rate_slow)); // site3->site2
    rates.push_back(pair<int,double>(4,rate_slow)); // site3->site4
    rates.push_back(pair<int,double>(7,rate_fast)); // site3->site7

    rates.push_back(pair<int,double>(3,rate_slow)); // site4->site3
    rates.push_back(pair<int,double>(8,rate_slow)); // site4->site8
    
    rates.push_back(pair<int,double>(1,rate_slow)); // site5->site1
    rates.push_back(pair<int,double>(6,rate_fast)); // site5->site6
    rates.push_back(pair<int,double>(9,rate_slow)); // site5->site9

    rates.push_back(pair<int,double>(2,rate_slow)); // site6->site2
    rates.push_back(pair<int,double>(5,rate_slow)); // site6->site5
    rates.push_back(pair<int,double>(7,rate_fast)); // site6->site7 
    rates.push_back(pair<int,double>(10,rate_slow)); // site6->site10

    rates.push_back(pair<int,double>(3,rate_slow)); // site7->site3
    rates.push_back(pair<int,double>(6,rate_fast)); // site7->site6
    rates.push_back(pair<int,double>(8,rate_slow)); // site7->site8
    rates.push_back(pair<int,double>(11,rate_slow)); // site7->site11

    rates.push_back(pair<int,double>(4,rate_slow)); // site8->site4
    rates.push_back(pair<int,double>(7,rate_fast)); // site8->site7
    rates.push_back(pair<int,double>(12,rate_slow)); // site8->site12

    rates.push_back(pair<int,double>(5,rate_slow)); // site9->site5
    rates.push_back(pair<int,double>(10,rate_slow)); // site9->site10

    rates.push_back(pair<int,double>(6,rate_fast)); // site10->site6
    rates.push_back(pair<int,double>(9,rate_slow)); // site10->site9
    rates.push_back(pair<int,double>(11,rate_slow)); // site10->site11

    rates.push_back(pair<int,double>(7,rate_fast)); // site11->site7
    rates.push_back(pair<int,double>(10,rate_slow)); // site11->site10
    rates.push_back(pair<int,double>(12,rate_slow)); // site11->site12

    rates.push_back(pair<int,double>(8,rate_slow)); // site12->site8
    rates.push_back(pair<int,double>(11,rate_slow)); // site12->site11

    // Ids of each of the sites
    vector<int> idsOfEachSite = { 
      1,  2,  3,  4,
      5,  6,  7,  8,
      9, 10, 11, 12};

    // Number of neighbors each site has
    vector<int> numberOfNeighbors = {
      2, 3, 3, 2,
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

    for( int index=0;index<12;++index){

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

    CoarseGrainSystem CGsystem;
    CGsystem.setTimeResolution(10.0);
    CGsystem.initializeSystem(ratesToNeighbors);
  }

  cout << "Testing: hop" << endl;
  {
    // In this example we will define 12 sites
    //
    // site1 - site2 - site3 - site4
    //   |       |       |       |
    // site5 - site6 - site7 - site8
    //   |       |       |       |
    // site9 - site10- site11- site12
    
    // Define a vector with all the rates starting from site 1 and describing
    // the rates to its neighbors and then continue with the other sites
    
    double rate_fast = 100;
    double rate_slow = 1;
    double rate_very_slow = 0.001;

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

    rates.push_back(pair<int,double>(2,rate_very_slow)); // site6->site2
    rates.push_back(pair<int,double>(5,rate_very_slow)); // site6->site5
    rates.push_back(pair<int,double>(7,rate_fast)); // site6->site7 
    rates.push_back(pair<int,double>(10,rate_very_slow)); // site6->site10

    rates.push_back(pair<int,double>(3,rate_very_slow)); // site7->site3
    rates.push_back(pair<int,double>(6,rate_fast)); // site7->site6
    rates.push_back(pair<int,double>(8,rate_very_slow)); // site7->site8
    rates.push_back(pair<int,double>(11,rate_very_slow)); // site7->site11

    rates.push_back(pair<int,double>(4,rate_slow)); // site8->site4
    rates.push_back(pair<int,double>(7,rate_slow)); // site8->site7
    rates.push_back(pair<int,double>(12,rate_slow)); // site8->site12

    rates.push_back(pair<int,double>(5,rate_slow)); // site9->site5
    rates.push_back(pair<int,double>(10,rate_slow)); // site9->site10

    rates.push_back(pair<int,double>(6,rate_slow)); // site10->site6
    rates.push_back(pair<int,double>(9,rate_slow)); // site10->site9
    rates.push_back(pair<int,double>(11,rate_slow)); // site10->site11

    rates.push_back(pair<int,double>(7,rate_slow)); // site11->site7
    rates.push_back(pair<int,double>(10,rate_slow)); // site11->site10
    rates.push_back(pair<int,double>(12,rate_slow)); // site11->site12

    rates.push_back(pair<int,double>(8,rate_slow)); // site12->site8
    rates.push_back(pair<int,double>(11,rate_slow)); // site12->site11

    // Ids of each of the sites
    vector<int> idsOfEachSite = { 
      1,  2,  3,  4,
      5,  6,  7,  8,
      9, 10, 11, 12};

    // Number of neighbors each site has
    vector<int> numberOfNeighbors = {
      2, 3, 3, 2,
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

    for( int index=0;index<12;++index){

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

    double time_limit = 10000;
    cout << "Running without cluster" << endl;
    // Without cluster
    {
      CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(1);
      double time_resolution = time_limit/10.0;
      CGsystem.setTimeResolution(time_resolution);
      CGsystem.setMinCoarseGrainIterationThreshold(constants::inf_iterations);
      CGsystem.initializeSystem(ratesToNeighbors);

      class Electron : public Walker {};

      Electron electron;
      // Place the electron on site 1
      int siteId = 1;
      electron.occupySite(siteId);

      vector<pair<int,Walker>> electrons;
      electrons.push_back(pair<int,Walker>(1,electron));

      CGsystem.initializeWalkers(electrons);

      vector<double> time_spent_on_sites(12,0.0);
      vector<int> hops_made_to_sites(12,0);

      double time = 0.0;
      int hop_count = 0;
      Walker& electron1 = electrons.at(0).second;
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
      CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(1);
      double time_resolution = time_limit/10.0;
      CGsystem.setTimeResolution(time_resolution);
      CGsystem.setPerformanceRatio(1.0);
      CGsystem.setMinCoarseGrainIterationThreshold(1000);
      CGsystem.initializeSystem(ratesToNeighbors);

      class Electron : public Walker {};

      Electron electron;
      // Place the electron on site 1
      int siteId = 1;
      electron.occupySite(siteId);

      vector<pair<int,Walker>> electrons;
      electrons.push_back(pair<int,Walker>(0,electron));

      CGsystem.initializeWalkers(electrons);

      vector<double> time_spent_on_sites(12,0.0);
      vector<int> hops_made_to_sites(12,0);

      double time = 0.0;
      int hop_count = 0;
      Walker& electron1 = electrons.at(0).second;
      int id = electrons.at(0).first;
      while(time<time_limit){
        CGsystem.hop(id,electron1);
        time_spent_on_sites.at(electron1.getIdOfSiteCurrentlyOccupying()-1) =
          electron1.getDwellTime();
        hops_made_to_sites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)++;
        time += electron1.getDwellTime();
        ++hop_count;
      }    

      auto clusters = CGsystem.getClusters();
      assert(clusters.size()==1);
      bool site6_found = false;
      bool site7_found = false;
      for( auto siteId : clusters.at(0) ){
        if(siteId==6) site6_found = true;
        if(siteId==7) site7_found = true;
      }
      assert(site6_found);
      assert(site7_found);
      // Check that the appropriate sites have been found in the cluster

      double sum_times = 0.0;
      for(auto time_site : time_spent_on_sites) sum_times+=time_site;

      double sum_time_ratio = 0.0;
      vector<double> time_ratios;
      for(auto time_site : time_spent_on_sites){
        time_ratios.push_back(time_site/sum_times);
        sum_time_ratio +=time_site/sum_times;
      }

    } // With cluster formation
  }

  cout << "Testing: hop 2" << endl;
  {
    // In this example we will define 14 sites
    //
    //         neigh6 - neigh7- neigh8
    //            |       |       |
    // neigh14- site1 - site2 - site3 - neigh9
    //            |       |       |
    // neigh13- site5 - site4 - neigh10
    //            |       |        
    //         neigh12- neigh11        
    
    // Define a vector with all the rates starting from site 1 and describing
    // the rates to its neighbors and then continue with the other sites
    //
    // We will exclude rates going from the neighbors to the inner sites
    
    double rate1 = 1.0/1.0;
    double rate2 = 1.0/400.0;
    double rate3 = 1.0/100.0;
    double rate4 = 1.0/25.0;
    double rate5 = 1.0/300.0;
    double rate6 = 1.0/0.1;

    vector<pair<int,double>> rates;
    rates.push_back(pair<int,double>(14,rate2)); // site1->neigh14
    rates.push_back(pair<int,double>(6,rate3)); // site1->neigh6
    rates.push_back(pair<int,double>(2,rate1)); // site1->site2
    rates.push_back(pair<int,double>(5,rate1)); // site1->site5

    rates.push_back(pair<int,double>(1,rate1)); // site2->site1
    rates.push_back(pair<int,double>(7,rate2)); // site2->neigh7
    rates.push_back(pair<int,double>(3,rate1)); // site2->site3
    rates.push_back(pair<int,double>(4,rate6)); // site2->site4

    rates.push_back(pair<int,double>(2,rate4)); // site3->site2
    rates.push_back(pair<int,double>(8,rate2)); // site3->neigh8
    rates.push_back(pair<int,double>(9,rate2)); // site3->neigh9
    rates.push_back(pair<int,double>(10,rate2)); // site3->neigh10

    rates.push_back(pair<int,double>(5,rate1)); // site4->site5
    rates.push_back(pair<int,double>(2,rate1)); // site4->site2
    rates.push_back(pair<int,double>(10,rate5)); // site4->neigh10
    rates.push_back(pair<int,double>(11,rate2)); // site4->neigh11

    rates.push_back(pair<int,double>(13,rate3)); // site5->neigh13
    rates.push_back(pair<int,double>(1,rate1)); // site5->site1
    rates.push_back(pair<int,double>(4,rate1)); // site5->site4
    rates.push_back(pair<int,double>(12,rate2)); // site5->neigh12

    rates.push_back(pair<int,double>(1,rate2));
    rates.push_back(pair<int,double>(2,rate2));
    rates.push_back(pair<int,double>(3,rate2));
    rates.push_back(pair<int,double>(3,rate2));
    rates.push_back(pair<int,double>(3,rate2));
    rates.push_back(pair<int,double>(4,rate2));
    rates.push_back(pair<int,double>(4,rate2));
    rates.push_back(pair<int,double>(5,rate2));
    rates.push_back(pair<int,double>(5,rate2));
    rates.push_back(pair<int,double>(1,rate2));

    // Ids of each of the sites
    vector<int> idsOfEachSite = { 
      1,  2,  3,  4,
      5,  6,  7,  8,
      9, 10, 11, 12,
     13, 14};

    // Number of neighbors each site has
    vector<int> numberOfNeighbors = {
      4, 4, 4, 4,
      4, 1, 1, 1,
      1, 1, 1, 1,
      1, 1};

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

    for( int index=0;index<14;++index){

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

    cout << endl;
    double time_limit = 10000;
    cout << "Without Coarse graining" << endl;
   
    // Number of electrons used for both the following crude and coarse grained
    // simulation runs
    int NumberElectrons = 8000;
    int number_of_sites = 14;
    // Will be compared with Coarse grained version
    vector<double> probabilityOnNeighCrude; 
    // Without cluster formation
    vector<double> portionOfTimeOnSiteNoCluster(5,0.0);
    vector<double> hops_to_sites_no_cluster(number_of_sites,0);  
    {
      CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(1);
      double time_resolution = time_limit/10.0;
      CGsystem.setTimeResolution(time_resolution);
      CGsystem.setMinCoarseGrainIterationThreshold(constants::inf_iterations);
      CGsystem.initializeSystem(ratesToNeighbors);

      class Electron : public Walker {};

      // Store the number of hops to each site 1-14
      vector<double> timeOnSites(number_of_sites,0.0);
      // Store the escape time from the cluster for each electron
      vector<double> escapeTimes;


      for(int i=0; i<NumberElectrons;++i){
        Electron electron;
        int id = 1;
        // Alternate placing electrons on sites 1-5
        int initialSite =  (i%5)+1;
        electron.occupySite(initialSite);
        vector<pair<int,Walker>> electrons;
        electrons.push_back(pair<int,Walker>(id,electron));
        CGsystem.initializeWalkers(electrons);
        double totalTimeOnCluster = 0.0;
        Walker & electron1 = electrons.at(0).second;

        while(electron1.getIdOfSiteCurrentlyOccupying()<6){
          CGsystem.hop(id,electron1);
          timeOnSites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)+=electron1.getDwellTime();
          totalTimeOnCluster+=electron1.getDwellTime(); 
        }

        CGsystem.removeWalkerFromSystem(id,electron1);
        escapeTimes.push_back(totalTimeOnCluster);
      }

      cout << "Total number of visits to each site" << endl;
      for(int site_id = 1; site_id <= number_of_sites; ++site_id){
        int visits = CGsystem.getVisitFrequencyOfSite(site_id);
        hops_to_sites_no_cluster.at(site_id-1) = static_cast<double>(visits); 
        cout << "id: " << site_id << " visits " << static_cast<double>(visits) << endl;
      }
    
      auto clusters = CGsystem.getClusters();
      // There should be no clusters found because the Threshold is so high
      assert(clusters.size()==0);

      int totalneigh = 0;
      totalneigh+= hops_to_sites_no_cluster.at(5);
      totalneigh+= hops_to_sites_no_cluster.at(6);
      totalneigh+= hops_to_sites_no_cluster.at(7);
      totalneigh+= hops_to_sites_no_cluster.at(8);
      totalneigh+= hops_to_sites_no_cluster.at(9);
      totalneigh+= hops_to_sites_no_cluster.at(10);
      totalneigh+= hops_to_sites_no_cluster.at(11);
      totalneigh+= hops_to_sites_no_cluster.at(12);
      totalneigh+= hops_to_sites_no_cluster.at(13);

      double value = static_cast<double>(hops_to_sites_no_cluster.at(5))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(6))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(7))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(8))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(9))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(10))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(11))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(12))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);
      value = static_cast<double>(hops_to_sites_no_cluster.at(13))/static_cast<double>(totalneigh);
      probabilityOnNeighCrude.push_back(value);

      cout << "Probability of Hopping to a neighboring site should match" << endl;
      for(int i=0; i<9;++i){
        cout << "Probability hop to neigh " << (i+6) << " " << probabilityOnNeighCrude.at(i) << endl;
      }

      assert(probabilityOnNeighCrude.at(0)<0.11);
      assert(probabilityOnNeighCrude.at(0)>0.09);

      assert(probabilityOnNeighCrude.at(1)<0.010);
      assert(probabilityOnNeighCrude.at(1)>0.003);

      assert(probabilityOnNeighCrude.at(2)<0.19);
      assert(probabilityOnNeighCrude.at(2)>0.15);

      assert(probabilityOnNeighCrude.at(3)<0.185);
      assert(probabilityOnNeighCrude.at(3)>0.16);

      assert(probabilityOnNeighCrude.at(4)<0.265);
      assert(probabilityOnNeighCrude.at(4)>0.23);

      assert(probabilityOnNeighCrude.at(5)<0.065);
      assert(probabilityOnNeighCrude.at(5)>0.05);

      assert(probabilityOnNeighCrude.at(6)<0.05);
      assert(probabilityOnNeighCrude.at(6)>0.03);

      assert(probabilityOnNeighCrude.at(7)<0.18);
      assert(probabilityOnNeighCrude.at(7)>0.16);

      assert(probabilityOnNeighCrude.at(8)<0.03);
      assert(probabilityOnNeighCrude.at(8)>0.02);

      double totalTimeOnSites = 0.0;
      totalTimeOnSites+=timeOnSites.at(0);
      totalTimeOnSites+=timeOnSites.at(1);
      totalTimeOnSites+=timeOnSites.at(2);
      totalTimeOnSites+=timeOnSites.at(3);
      totalTimeOnSites+=timeOnSites.at(4);

      portionOfTimeOnSiteNoCluster.at(0)=timeOnSites.at(0)/totalTimeOnSites;
      portionOfTimeOnSiteNoCluster.at(1)=timeOnSites.at(1)/totalTimeOnSites;
      portionOfTimeOnSiteNoCluster.at(2)=timeOnSites.at(2)/totalTimeOnSites;
      portionOfTimeOnSiteNoCluster.at(3)=timeOnSites.at(3)/totalTimeOnSites;
      portionOfTimeOnSiteNoCluster.at(4)=timeOnSites.at(4)/totalTimeOnSites;

      cout << "Time percentage on sites" << endl; 
      for(int i=0; i<5;++i){
        cout << "Percent time spent on site " << (i+1) << " ";
        cout << portionOfTimeOnSiteNoCluster.at(i) << endl;
      }

      assert(
          portionOfTimeOnSiteNoCluster.at(0)>0.075 && 
          portionOfTimeOnSiteNoCluster.at(0) < 0.09);
      assert(
          portionOfTimeOnSiteNoCluster.at(1)>0.02 && 
          portionOfTimeOnSiteNoCluster.at(1) < 0.03);
      assert(
          portionOfTimeOnSiteNoCluster.at(2)>0.54 && 
          portionOfTimeOnSiteNoCluster.at(2) < 0.56);
      assert(
          portionOfTimeOnSiteNoCluster.at(3)>0.19 && 
          portionOfTimeOnSiteNoCluster.at(3) < 0.21);
      assert(
          portionOfTimeOnSiteNoCluster.at(4)>0.13 && 
          portionOfTimeOnSiteNoCluster.at(4) < 0.15);
    } // Without Cluster formation

    cout << endl;
    cout << "With Coarse graining" << endl;
    // With cluster formation
    {

      // Store the number of hops to each site 1-14
      vector<double> timeOnSites(number_of_sites,0.0);
      // Store the escape time from the cluster for each electron
      vector<double> escapeTimes;

      CoarseGrainSystem CGsystem;
      CGsystem.setMinCoarseGrainIterationThreshold(10);
      CGsystem.setPerformanceRatio(0.2);
      CGsystem.setRandomSeed(1);
      double time_resolution = time_limit/10.0;
      CGsystem.setTimeResolution(time_resolution);
      CGsystem.initializeSystem(ratesToNeighbors);
      int cycles = 1;
      for(int cycle = 0; cycle < cycles ;++cycle ){
        class Electron : public Walker {};


        for(int i=0; i<NumberElectrons;++i){
          Electron electron;
          int id = 1;
          // Alternate placing electrons on sites 1-5
          int initialSite =  (i%5)+1;
          electron.occupySite(initialSite);
          vector<pair<int,Walker>> electrons;
          electrons.push_back(pair<int,Walker>(id,electron));
          CGsystem.initializeWalkers(electrons);
          double totalTimeOnCluster = 0.0;
          // First hop is ignored 
          Walker & electron1 = electrons.at(0).second;
          CGsystem.hop(id,electron1);

          while(electron1.getIdOfSiteCurrentlyOccupying()<6){
            timeOnSites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)+=electron1.getDwellTime();
            totalTimeOnCluster+=electron1.getDwellTime(); 
            CGsystem.hop(id,electron1);
          }

          CGsystem.removeWalkerFromSystem(id,electron1);
          escapeTimes.push_back(totalTimeOnCluster);
        }

      }

      cout << "Total number of visits to each site" << endl;
      vector<double> hops_to_sites(number_of_sites,0);  
      for(int site_id = 1; site_id <= number_of_sites; ++site_id){
        int visits = CGsystem.getVisitFrequencyOfSite(site_id);
        hops_to_sites.at(site_id-1) = static_cast<double>(visits)/static_cast<double>(cycles); 
        cout << "id: " << site_id << " visits " << static_cast<double>(visits)/static_cast<double>(cycles) << endl;
      }

      assert(hops_to_sites.at(0)<hops_to_sites_no_cluster.at(0)*1.2);
      assert(hops_to_sites.at(0)>hops_to_sites_no_cluster.at(0)*0.8);
      assert(hops_to_sites.at(1)<hops_to_sites_no_cluster.at(1)*1.2);
      assert(hops_to_sites.at(1)>hops_to_sites_no_cluster.at(1)*0.8);
      assert(hops_to_sites.at(2)<hops_to_sites_no_cluster.at(2)*1.2);
      assert(hops_to_sites.at(2)>hops_to_sites_no_cluster.at(2)*0.8);
      assert(hops_to_sites.at(3)<hops_to_sites_no_cluster.at(3)*1.2);
      assert(hops_to_sites.at(3)>hops_to_sites_no_cluster.at(3)*0.8);
      assert(hops_to_sites.at(4)<hops_to_sites_no_cluster.at(4)*1.2);
      assert(hops_to_sites.at(4)>hops_to_sites_no_cluster.at(4)*0.8);
      assert(hops_to_sites.at(5)<hops_to_sites_no_cluster.at(5)*1.2);
      assert(hops_to_sites.at(5)>hops_to_sites_no_cluster.at(5)*0.8);
      assert(hops_to_sites.at(6)<hops_to_sites_no_cluster.at(6)*1.2);
      assert(hops_to_sites.at(6)>hops_to_sites_no_cluster.at(6)*0.8);
      assert(hops_to_sites.at(7)<hops_to_sites_no_cluster.at(7)*1.2);
      assert(hops_to_sites.at(7)>hops_to_sites_no_cluster.at(7)*0.8);
      assert(hops_to_sites.at(8)<hops_to_sites_no_cluster.at(8)*1.2);
      assert(hops_to_sites.at(8)>hops_to_sites_no_cluster.at(8)*0.8);
      assert(hops_to_sites.at(9)<hops_to_sites_no_cluster.at(9)*1.2);
      assert(hops_to_sites.at(9)>hops_to_sites_no_cluster.at(9)*0.8);
      assert(hops_to_sites.at(10)<hops_to_sites_no_cluster.at(10)*1.2);
      assert(hops_to_sites.at(10)>hops_to_sites_no_cluster.at(10)*0.8);
      assert(hops_to_sites.at(11)<hops_to_sites_no_cluster.at(11)*1.2);
      assert(hops_to_sites.at(11)>hops_to_sites_no_cluster.at(11)*0.8);
      assert(hops_to_sites.at(12)<hops_to_sites_no_cluster.at(12)*1.2);
      assert(hops_to_sites.at(12)>hops_to_sites_no_cluster.at(12)*0.8);
      assert(hops_to_sites.at(13)<hops_to_sites_no_cluster.at(13)*1.2);
      assert(hops_to_sites.at(13)>hops_to_sites_no_cluster.at(13)*0.8);

      auto clusters = CGsystem.getClusters();
      cout << "Clusters size " << clusters.size() << endl;
      assert(clusters.size()==1);
      bool site1_found = false;
      bool site2_found = false;
      bool site3_found = false;
      bool site4_found = false;
      bool site5_found = false;
      for( auto siteId : clusters.begin()->second ){
        if(siteId==1) site1_found = true;
        if(siteId==2) site2_found = true;
        if(siteId==3) site3_found = true;
        if(siteId==4) site4_found = true;
        if(siteId==5) site5_found = true;
      }
      assert(site1_found);
      assert(site2_found);
      assert(site3_found);
      assert(site4_found);
      assert(site5_found);


      int totalneigh = 0;
      totalneigh+= hops_to_sites.at(5);
      totalneigh+= hops_to_sites.at(6);
      totalneigh+= hops_to_sites.at(7);
      totalneigh+= hops_to_sites.at(8);
      totalneigh+= hops_to_sites.at(9);
      totalneigh+= hops_to_sites.at(10);
      totalneigh+= hops_to_sites.at(11);
      totalneigh+= hops_to_sites.at(12);
      totalneigh+= hops_to_sites.at(13);

      vector<double> probabilityOnNeigh; 
      double value = static_cast<double>(hops_to_sites.at(5))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(6))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(7))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(8))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(9))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(10))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(11))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(12))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);
      value = static_cast<double>(hops_to_sites.at(13))/static_cast<double>(totalneigh);
      probabilityOnNeigh.push_back(value);

      cout << "Probability of Hopping to a neighboring site should match" << endl;
      for(int i=0; i<9;++i){
        cout << "Probability hop to neigh " << (i+6) << " " << probabilityOnNeigh.at(i) << endl;
      }

      assert(probabilityOnNeigh.at(0)<probabilityOnNeighCrude.at(0)*1.2);
      assert(probabilityOnNeigh.at(0)>probabilityOnNeighCrude.at(0)*0.8);

      assert(probabilityOnNeigh.at(1)<probabilityOnNeighCrude.at(1)*1.2);
      assert(probabilityOnNeigh.at(1)>probabilityOnNeighCrude.at(1)*0.8);

      assert(probabilityOnNeigh.at(2)<probabilityOnNeighCrude.at(2)*1.2);
      assert(probabilityOnNeigh.at(2)>probabilityOnNeighCrude.at(2)*0.8);

      assert(probabilityOnNeigh.at(3)<probabilityOnNeighCrude.at(3)*1.2);
      assert(probabilityOnNeigh.at(3)>probabilityOnNeighCrude.at(3)*0.8);

      assert(probabilityOnNeigh.at(4)<probabilityOnNeighCrude.at(4)*1.2);
      assert(probabilityOnNeigh.at(4)>probabilityOnNeighCrude.at(4)*0.8);

      assert(probabilityOnNeigh.at(5)<probabilityOnNeighCrude.at(5)*1.2);
      assert(probabilityOnNeigh.at(5)>probabilityOnNeighCrude.at(5)*0.8);

      assert(probabilityOnNeigh.at(6)<probabilityOnNeighCrude.at(6)*1.2);
      assert(probabilityOnNeigh.at(6)>probabilityOnNeighCrude.at(6)*0.8);

      assert(probabilityOnNeigh.at(7)<probabilityOnNeighCrude.at(7)*1.2);
      assert(probabilityOnNeigh.at(7)>probabilityOnNeighCrude.at(7)*0.8);

      assert(probabilityOnNeigh.at(8)<probabilityOnNeighCrude.at(8)*1.2);
      assert(probabilityOnNeigh.at(8)>probabilityOnNeighCrude.at(8)*0.8);

      double totalTimeOnSites = 0.0;
      totalTimeOnSites+=timeOnSites.at(0);
      totalTimeOnSites+=timeOnSites.at(1);
      totalTimeOnSites+=timeOnSites.at(2);
      totalTimeOnSites+=timeOnSites.at(3);
      totalTimeOnSites+=timeOnSites.at(4);

      vector<double> portionOfTimeOnSite(5,0.0);
      portionOfTimeOnSite.at(0)=timeOnSites.at(0)/totalTimeOnSites;
      portionOfTimeOnSite.at(1)=timeOnSites.at(1)/totalTimeOnSites;
      portionOfTimeOnSite.at(2)=timeOnSites.at(2)/totalTimeOnSites;
      portionOfTimeOnSite.at(3)=timeOnSites.at(3)/totalTimeOnSites;
      portionOfTimeOnSite.at(4)=timeOnSites.at(4)/totalTimeOnSites;

      cout << "Time percentage on sites" << endl; 
      for(int i=0; i<5;++i){
        cout << "Percent time spent on site " << (i+1) << " " << portionOfTimeOnSite.at(i) << endl;
      }

      assert(
          portionOfTimeOnSite.at(0)>portionOfTimeOnSiteNoCluster.at(0)*0.8 && 
          portionOfTimeOnSite.at(0)<portionOfTimeOnSiteNoCluster.at(0)*1.2);
      assert(
          portionOfTimeOnSite.at(1)>portionOfTimeOnSiteNoCluster.at(1)*0.8 && 
          portionOfTimeOnSite.at(1)<portionOfTimeOnSiteNoCluster.at(1)*1.2);
      assert(
          portionOfTimeOnSite.at(2)>portionOfTimeOnSiteNoCluster.at(2)*0.8 && 
          portionOfTimeOnSite.at(2)<portionOfTimeOnSiteNoCluster.at(2)*1.2);
      assert(
          portionOfTimeOnSite.at(3)>portionOfTimeOnSiteNoCluster.at(3)*0.8 && 
          portionOfTimeOnSite.at(3)<portionOfTimeOnSiteNoCluster.at(3)*1.2);
      assert(
          portionOfTimeOnSite.at(4)>portionOfTimeOnSiteNoCluster.at(4)*0.8 && 
          portionOfTimeOnSite.at(4)<portionOfTimeOnSiteNoCluster.at(4)*1.2);

    }// With cluster formation
  }

	return 0;
}
