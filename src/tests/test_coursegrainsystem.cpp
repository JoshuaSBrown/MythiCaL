#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "../../include/kmccoursegrain/coursegrainsystem.hpp"

using namespace std;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: CourseGrainSystem constructor" << endl;
  {
    CourseGrainSystem CGsystem;
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

    CourseGrainSystem CGsystem;
    CGsystem.initializeSystem(ratesToNeighbors);
  }
	return 0;
}
