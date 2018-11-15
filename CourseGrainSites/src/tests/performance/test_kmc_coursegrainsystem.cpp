#include <iostream>
#include <cassert>
#include <vector>
#include <memory>
#include <chrono>

#include "../../../include/kmccoursegrain/kmc_coursegrainsystem.hpp"
#include "../../../include/kmccoursegrain/kmc_particle.hpp"

using namespace std;
using namespace std::chrono;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: kmc_coursegraingsystem" << endl;
  cout << "This executable tests the coursegraining performance of " << endl;
  cout << "code when it is able to identify a cluster during the " << endl;
  cout << "runtime and compares it to the code when a cluster is not " << endl;
  cout << "identified." << endl;

 
  high_resolution_clock::time_point clusterStart; 
  high_resolution_clock::time_point clusterEnd; 
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
    unordered_map< int,unordered_map< int,double *>> ratesToNeighbors;

    int global_index = 0;

    for( int index=0;index<14;++index){

      unordered_map< int,double*> ratesFromSiteToNeighbors;

      for( int rate_index = 0; 
          rate_index<numberOfNeighbors.at(index);
          ++rate_index){

        ratesFromSiteToNeighbors[rates.at(global_index).first] = \
           &(rates.at(global_index).second);

        ++global_index;        
      }
  
      ratesToNeighbors[idsOfEachSite.at(index)] = ratesFromSiteToNeighbors;
    }

    
    KMC_CourseGrainSystem CGsystem;
    CGsystem.setRandomSeed(1);
    CGsystem.setCourseGrainIterationThreshold(10000);
    CGsystem.initializeSystem(ratesToNeighbors);
    
    class Electron : public KMC_Particle {};
    
    // Store the number of hops to each site 1-14
    vector<int> hopsToSites(14,0);  
    vector<double> timeOnSites(14,0.0);
    // Store the escape time from the cluster for each electron
    vector<double> escapeTimes;

    int NumberElectrons = 5000;
    cout << "Cluster performance test starting" << endl;
    clusterStart = high_resolution_clock::now();
    for(int i=0; i<NumberElectrons;++i){
      Electron electron;
      // Alternate placing electrons on sites 1-5
      int initialSite =  (i%5)+1;
      electron.occupySite(initialSite);
      vector<KMC_Particle> electrons;
      electrons.push_back(electron);
      CGsystem.initializeParticles(electrons);
      
      KMC_Particle & electron1 = electrons.at(0);
      double totalTimeOnCluster = 0.0;
      while(electron1.getIdOfSiteCurrentlyOccupying()<6){

        CGsystem.hop(electron1);
        hopsToSites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)++;
        timeOnSites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)+=electron1.getDwellTime();
        totalTimeOnCluster+=electron1.getDwellTime(); 
      }

      CGsystem.removeParticleFromSystem(electron1);
      escapeTimes.push_back(totalTimeOnCluster);
    }
    clusterEnd = high_resolution_clock::now();
  }

  high_resolution_clock::time_point siteStart; 
  high_resolution_clock::time_point siteEnd; 
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
    unordered_map< int,unordered_map< int,double *>> ratesToNeighbors;

    int global_index = 0;

    for( int index=0;index<14;++index){

      unordered_map< int,double*> ratesFromSiteToNeighbors;

      for( int rate_index = 0; 
          rate_index<numberOfNeighbors.at(index);
          ++rate_index){

        ratesFromSiteToNeighbors[rates.at(global_index).first] = \
           &(rates.at(global_index).second);

        ++global_index;        
      }
  
      ratesToNeighbors[idsOfEachSite.at(index)] = ratesFromSiteToNeighbors;
    }

    KMC_CourseGrainSystem CGsystem;
    CGsystem.setRandomSeed(1);
    CGsystem.setCourseGrainIterationThreshold(1000000);
    CGsystem.initializeSystem(ratesToNeighbors);
    
    class Electron : public KMC_Particle {};
    
    // Store the number of hops to each site 1-14
    vector<int> hopsToSites(14,0);  
    vector<double> timeOnSites(14,0.0);
    // Store the escape time from the cluster for each electron
    vector<double> escapeTimes;

    int NumberElectrons = 5000;

    cout << "Site performance test starting" << endl;
    siteStart = high_resolution_clock::now();
    for(int i=0; i<NumberElectrons;++i){
      Electron electron;
      // Alternate placing electrons on sites 1-5
      int initialSite =  (i%5)+1;
      electron.occupySite(initialSite);
      vector<KMC_Particle> electrons;
      electrons.push_back(electron);
      CGsystem.initializeParticles(electrons);
      double totalTimeOnCluster = 0.0;
      KMC_Particle & electron1 = electrons.at(0);
      while(electron1.getIdOfSiteCurrentlyOccupying()<6){
        CGsystem.hop(electron1);
        hopsToSites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)++;
        timeOnSites.at(electron1.getIdOfSiteCurrentlyOccupying()-1)+=electron1.getDwellTime();
        totalTimeOnCluster+=electron1.getDwellTime(); 
      }

      CGsystem.removeParticleFromSystem(electron1);
      escapeTimes.push_back(totalTimeOnCluster);
    }
    siteEnd = high_resolution_clock::now();
  }

  auto durationCluster = duration_cast<microseconds>(clusterEnd-clusterStart).count();
  auto durationSite = duration_cast<microseconds>(siteEnd-siteStart).count();

  cout << "Cluster time " << durationCluster << endl;
  cout << "Site time " << durationSite << endl;

  assert(durationCluster<durationSite);
  return 0;
}
