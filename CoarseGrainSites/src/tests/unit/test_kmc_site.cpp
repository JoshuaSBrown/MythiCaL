#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "../../libkmccoarsegrain/topologyfeatures/kmc_site.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

  cout << "Testing: Site constructor" << endl;
  {
    KMC_Site site;
  }

  cout << "Testing: Site Id setter" << endl;
  {
    KMC_Site site;
    site.setId(0);
  }
	
  cout << "Testing: Site Id getter" << endl;
  {
    KMC_Site site;
    site.setId(0);
    assert(site.getId()==0);

    bool fail = false;
    KMC_Site site2;
    try {
      site2.getId();
    }catch(...){
      fail = true;
    }
    assert(fail);
  }

  cout << "Testing: setRatesToNeighbors" << endl;
  {
    unordered_map< int, double> neighRates;
    double rate1 = 400;
    double rate2 = 200;
    double rate3 = 10;
    double rate4 = 1;
    neighRates[1]=rate1;
    neighRates[2]=rate2;
    neighRates[3]=rate3;
    neighRates[4]=rate4;
    
    KMC_Site site;
    site.setRatesToNeighbors(&neighRates);
  }

  cout << "Testing: getRateToNeighbor" << endl;
  {
    unordered_map< int, double> neighRates;
    double rate1 = 400;
    double rate2 = 200;
    double rate3 = 10;
    double rate4 = 1;
    neighRates[1]=rate1;
    neighRates[2]=rate2;
    neighRates[3]=rate3;
    neighRates[4]=rate4;
    
    KMC_Site site;
    site.setRatesToNeighbors(&neighRates);
    assert(static_cast<int>(site.getRateToNeighbor(1))==400);
    assert(static_cast<int>(site.getRateToNeighbor(2))==200);
    assert(static_cast<int>(site.getRateToNeighbor(3))==10);
    assert(static_cast<int>(site.getRateToNeighbor(4))==1);
  }
 /* 
  cout << "Testing: addNeighRate" << endl;
  {
    double rate1 = 400;
    double rate2 = 200;
    double rate3 = 10;
    double rate4 = 1;
    
    KMC_Site site;

    site.addNeighRate(pair< int,double *>(1,&rate1));
    site.addNeighRate(pair< int,double *>(2,&rate2));
    site.addNeighRate(pair< int,double *>(3,&rate3));
    site.addNeighRate(pair< int,double * >(4,&rate4));

    assert(static_cast<int>(site.getRateToNeighbor(1))==400);
    assert(static_cast<int>(site.getRateToNeighbor(2))==200);
    assert(static_cast<int>(site.getRateToNeighbor(3))==10);
    assert(static_cast<int>(site.getRateToNeighbor(4))==1);


  }*/

/*  cout << "Testing: resetNeighRate" << endl;
  {
    double rate1 = 400;
    double rate2 = 200;
    double rate3 = 10;
    double rate4 = 1;
    
    KMC_Site site;

    site.resetNeighRate(pair< int,double * >(1,&rate1));
    site.resetNeighRate(pair< int,double * >(2,&rate2));
    site.resetNeighRate(pair< int,double * >(3,&rate3));
    site.resetNeighRate(pair< int,double * >(4,&rate4));

    int value = static_cast<int>(site.getRateToNeighbor(1));
    assert(value==400);
    value = static_cast<int>(site.getRateToNeighbor(2));
    assert(value ==200);
    value = static_cast<int>(site.getRateToNeighbor(3));
    assert(value==10);
    value = static_cast<int>(site.getRateToNeighbor(4));
    assert(value==1);

    site.resetNeighRate(pair< int,double * >(1,&rate4));
    site.resetNeighRate(pair< int,double * >(2,&rate3));
    site.resetNeighRate(pair< int,double * >(3,&rate2));
    site.resetNeighRate(pair< int,double * >(4,&rate1));
 
    value = static_cast<int>(site.getRateToNeighbor(4));
    assert(value==400);
    value = static_cast<int>(site.getRateToNeighbor(3));
    assert(value==200);
    value = static_cast<int>(site.getRateToNeighbor(2));
    assert(value==10);
    value =static_cast<int>(site.getRateToNeighbor(1));
    assert(value==1);
  }*/

  cout << "Testing: isNeighbor" << endl;
  {
    double rate1 = 400;
    double rate2 = 200;
    double rate3 = 10;
    double rate4 = 1;
   
    int id = 0; 
    KMC_Site site(id);

    unordered_map<int,double> rates;
    rates[1] = rate1;
    rates[2] = rate2;
    rates[3] = rate3;
    rates[4] = rate4;

    site.setRatesToNeighbors(&rates);
/*    site.addNeighRate(pair< int,double * >(1,&rate1));
    site.addNeighRate(pair< int,double * >(2,&rate2));
    site.addNeighRate(pair< int,double * >(3,&rate3));
    site.addNeighRate(pair< int,double * >(4,&rate4));
*/
    assert(site.isNeighbor(0)==false);
    assert(site.isNeighbor(1)==true);
    assert(site.isNeighbor(2)==true);
    assert(site.isNeighbor(3)==true);
    assert(site.isNeighbor(4)==true);
    assert(site.isNeighbor(5)==false);
  }

  cout << "Testing: occupation functions" << endl;
  {
    KMC_Site site;
    assert(site.isOccupied()==false);
    site.occupy();
    assert(site.isOccupied()==true);
    site.vacate();
    assert(site.isOccupied()==false);
  }

  cout << "Testing: cluster functions " << endl;
  {
    KMC_Site site;
    assert(site.partOfCluster()==false);
    site.setClusterId(1);
    assert(site.partOfCluster()==true);
    assert(site.getClusterId()==1);
  }

  cout << "Testing: probability to hop to neighbor" << endl;
  {
    unordered_map< int, double > neighRates;
    double rate1 = 1;
    double rate2 = 1;
    double rate3 = 1;
    double rate4 = 1;
    neighRates[1]=rate1;
    neighRates[2]=rate2;
    neighRates[3]=rate3;
    neighRates[4]=rate4;
    
    KMC_Site site;
    site.setRatesToNeighbors(&neighRates);

/*    bool fail = false;
    try {
      site.getProbabilityOfHoppingToNeighboringSite(0);
    } catch(...) {
      fail = true;
    }
    assert(fail);
*/
    double probability = site.getProbabilityOfHoppingToNeighboringSite(1);
    cout << "probability " << probability << endl;
    assert(static_cast<int>(probability*100)==25);
  
  }

  cout << "Testing: getNeighborSiteIds" << endl;
  {
    unordered_map< int, double > neighRates;
    double rate1 = 1;
    double rate2 = 1;
    double rate3 = 1;
    double rate4 = 1;
    neighRates[1]=rate1;
    neighRates[2]=rate2;
    neighRates[3]=rate3;
    neighRates[4]=rate4;
    
    KMC_Site site;
    site.setRatesToNeighbors(&neighRates);
    
    auto neighborIds = site.getNeighborSiteIds();

    bool site1_identified = false;
    bool site2_identified = false;
    bool site3_identified = false;
    bool site4_identified = false;
    
    for(auto neighId : neighborIds){
      if(neighId==1) site1_identified=true;
      if(neighId==2) site2_identified=true;
      if(neighId==3) site3_identified=true;
      if(neighId==4) site4_identified=true;
    }

    assert(site1_identified);
    assert(site2_identified);
    assert(site3_identified);
    assert(site4_identified);
  }

  cout << "Testing: site output" << endl;
  {
    unordered_map< int, double > neighRates;
    double rate1 = 1;
    double rate2 = 1;
    double rate3 = 1;
    double rate4 = 1;
    neighRates[1]=rate1;
    neighRates[2]=rate2;
    neighRates[3]=rate3;
    neighRates[4]=rate4;
    
    KMC_Site site;
    site.setId(0);
    site.setRatesToNeighbors(&neighRates);
    cout << site << endl;

  }

	return 0;
}
