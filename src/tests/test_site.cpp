#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include <kmccoursegrain/site.hpp>

using namespace std;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: Site constructor" << endl;
  {
    Site site;
  }

  cout << "Testing: Site Id setter" << endl;
  {
    Site site;
    site.setId(0);
  }
	
  cout << "Testing: Site Id getter" << endl;
  {
    Site site;
    site.setId(0);
    assert(site.getId()==0);

    bool fail = false;
    Site site2;
    try {
      site2.getId();
    }catch(...){
      fail = true;
    }
    assert(fail);
  }

  cout << "Testing: setting rates" << endl;
  {
    map<const int, double *> neighRates;
    double rate1 = 400;
    double rate2 = 200;
    double rate3 = 10;
    double rate4 = 1;
    neighRates[1]=&rate1;
    neighRates[2]=&rate2;
    neighRates[3]=&rate3;
    neighRates[4]=&rate4;
    
    Site site;
    site.setRatesToNeighbors(neighRates);
  }

  cout << "Testing: probability to hop to neighbor" << endl;
  {
    map<const int, double *> neighRates;
    double rate1 = 1;
    double rate2 = 1;
    double rate3 = 1;
    double rate4 = 1;
    neighRates[1]=&rate1;
    neighRates[2]=&rate2;
    neighRates[3]=&rate3;
    neighRates[4]=&rate4;
    
    Site site;
    site.setRatesToNeighbors(neighRates);

    bool fail = false;
    try {
      site.probHopToNeigh(0);
    } catch(...) {
      fail = true;
    }
    assert(fail);

    double probability = site.probHopToNeigh(1);
    assert(static_cast<int>(probability*100)==25);
  
  }

  cout << "Testing: site output" << endl;
  {
    map<const int, double *> neighRates;
    double rate1 = 1;
    double rate2 = 1;
    double rate3 = 1;
    double rate4 = 1;
    neighRates[1]=&rate1;
    neighRates[2]=&rate2;
    neighRates[3]=&rate3;
    neighRates[4]=&rate4;
    
    Site site;
    site.setId(0);
    site.setRatesToNeighbors(neighRates);
    cout << site << endl;

  }

	return 0;
}
