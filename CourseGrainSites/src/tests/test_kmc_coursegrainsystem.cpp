#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "../../include/kmccoursegrain/kmc_coursegrainsystem.hpp"
#include "../../include/kmccoursegrain/kmc_particle.hpp"

using namespace std;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: CourseGrainSystem constructor" << endl;
  {
    KMC_CourseGrainSystem CGsystem;
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
    map<const int,map<const int,double *>> ratesToNeighbors;

    int global_index = 0;

    for( int index=0;index<12;++index){

      map<const int,double*> ratesFromSiteToNeighbors;

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
    map<const int,map<const int,double *>> ratesToNeighbors;

    int global_index = 0;

    for( int index=0;index<12;++index){

      map<const int,double*> ratesFromSiteToNeighbors;

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
    CGsystem.initializeSystem(ratesToNeighbors);

    class Electron : public KMC_Particle {};
  
    Electron electron;
    // Place the electron on site 1
    electron.occupySite(1);

    auto electron_ptr = make_shared<Electron>(electron);
    vector<shared_ptr<KMC_Particle>> electrons;
    electrons.push_back(electron_ptr);

    CGsystem.initializeParticles(electrons);

    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==1);
    assert(electron_ptr->getPotentialSite()==5);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==34);

    CGsystem.hop(electron_ptr);
    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==5);
    assert(electron_ptr->getPotentialSite()==6);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==1);

    CGsystem.hop(electron_ptr);
    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==6);
    assert(electron_ptr->getPotentialSite()==7);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==15);

    CGsystem.hop(electron_ptr);
    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==7);
    assert(electron_ptr->getPotentialSite()==6);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==11);

    // Enough to trigger the formation of a cluster
    for(int i=0; i<50;++i){
      CGsystem.hop(electron_ptr);
    }    

    cout << static_cast<int>(electron_ptr->getDwellTime()*1000) << endl;
    CGsystem.hop(electron_ptr);
    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==7);
    assert(electron_ptr->getPotentialSite()==7);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==183);

    CGsystem.hop(electron_ptr);
    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==7);
    assert(electron_ptr->getPotentialSite()==7);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==120);

    CGsystem.hop(electron_ptr);
    assert(electron_ptr->getIdOfSiteCurrentlyOccupying()==7);
    assert(electron_ptr->getPotentialSite()==6);
    assert(static_cast<int>(electron_ptr->getDwellTime()*1000)==53);
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
    map<const int,map<const int,double *>> ratesToNeighbors;

    int global_index = 0;

    for( int index=0;index<14;++index){

      map<const int,double*> ratesFromSiteToNeighbors;

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
    CGsystem.initializeSystem(ratesToNeighbors);

    class Electron : public KMC_Particle {};
    
    // Store the number of hops to each site 1-14
    vector<int> hopsToSites(14,0);  
    vector<double> timeOnSites(14,0.0);
    // Store the escape time from the cluster for each electron
    vector<double> escapeTimes;

    int NumberElectrons = 1000;

    for(int i=0; i<NumberElectrons;++i){
      Electron electron;
      electron.setMemoryCapacity(4);
      // Alternate placing electrons on sites 1-5
      int initialSite =  (i%5)+1;
      electron.occupySite(initialSite);
      auto electron_ptr = make_shared<Electron>(electron);
      vector<shared_ptr<KMC_Particle>> electrons;
      electrons.push_back(electron_ptr);
      CGsystem.initializeParticles(electrons);
      double totalTimeOnCluster = 0.0;
      while(electron_ptr->getIdOfSiteCurrentlyOccupying()<6){
        CGsystem.hop(electron_ptr);
        hopsToSites.at(electron_ptr->getIdOfSiteCurrentlyOccupying()-1)++;
        timeOnSites.at(electron_ptr->getIdOfSiteCurrentlyOccupying()-1)+=electron_ptr->getDwellTime();
        totalTimeOnCluster+=electron_ptr->getDwellTime(); 
      }
      CGsystem.removeParticleFromSystem(electron_ptr);
      escapeTimes.push_back(totalTimeOnCluster);
    }
   
    int totalsites = 0;
    totalsites+= hopsToSites.at(0);
    totalsites+= hopsToSites.at(1);
    totalsites+= hopsToSites.at(2);
    totalsites+= hopsToSites.at(3);
    totalsites+= hopsToSites.at(4);

    vector<double> probabilityOnSite; 

    double value = static_cast<double>(hopsToSites.at(0))/static_cast<double>(totalsites);
    probabilityOnSite.push_back(value);
    value = static_cast<double>(hopsToSites.at(1))/static_cast<double>(totalsites);
    probabilityOnSite.push_back(value);
    value = static_cast<double>(hopsToSites.at(2))/static_cast<double>(totalsites);
    probabilityOnSite.push_back(value);
    value = static_cast<double>(hopsToSites.at(3))/static_cast<double>(totalsites);
    probabilityOnSite.push_back(value);
    value = static_cast<double>(hopsToSites.at(4))/static_cast<double>(totalsites);
    probabilityOnSite.push_back(value);

    for(int i=0; i<5;++i){
      cout << "Probability hop to site " << (i+1) << " " << probabilityOnSite.at(i) << endl;
    }
    assert(probabilityOnSite.at(0)>0.139);
    assert(probabilityOnSite.at(0)<0.142);
   
    assert(probabilityOnSite.at(1)>0.258); 
    assert(probabilityOnSite.at(1)<0.262); 

    assert(probabilityOnSite.at(2)>0.020);
    assert(probabilityOnSite.at(2)<0.023);

    assert(probabilityOnSite.at(3)>0.336);
    assert(probabilityOnSite.at(3)<0.338);

    assert(probabilityOnSite.at(4)>0.237);
    assert(probabilityOnSite.at(4)<0.251);

    int totalneigh = 0;
    totalneigh+= hopsToSites.at(5);
    totalneigh+= hopsToSites.at(6);
    totalneigh+= hopsToSites.at(7);
    totalneigh+= hopsToSites.at(8);
    totalneigh+= hopsToSites.at(9);
    totalneigh+= hopsToSites.at(10);
    totalneigh+= hopsToSites.at(11);
    totalneigh+= hopsToSites.at(12);
    totalneigh+= hopsToSites.at(13);

    vector<double> probabilityOnNeigh; 
    value = static_cast<double>(hopsToSites.at(5))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(6))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(7))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(8))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(9))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(10))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(11))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(12))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);
    value = static_cast<double>(hopsToSites.at(13))/static_cast<double>(totalneigh);
    probabilityOnNeigh.push_back(value);

    cout << "Probability of Hopping to a neighboring site should match" << endl;
    for(int i=0; i<9;++i){
      cout << "Probability hop to neigh " << (i+6) << " " << probabilityOnNeigh.at(i) << endl;
    }


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
  }

	return 0;
}
