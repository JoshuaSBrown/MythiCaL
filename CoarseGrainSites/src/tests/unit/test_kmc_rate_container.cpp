
#include <iostream>
#include <cassert>

#include "../../libkmccoarsegrain/kmc_rate_container.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){
  cout << "Testing: Constructor" << endl;
  {
    KMC_Rate_Container rate_container;
    Rate_Map rate_map;
    double rate = 1.0;
    rate_map[1][2] = &rate; 
    KMC_Rate_Container rate_containter2(rate_map);
  }

  cout << "Testing: incomingRateCount" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;

    KMC_Rate_Container rate_container(rate_map);
    assert(rate_container.incomingRateCount(2)==2);
    assert(rate_container.incomingRateCount(1)==0);
    assert(rate_container.incomingRateCount(3)==0);
  }
 
  cout << "Testing: outgoingRateCount" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;

    KMC_Rate_Container rate_container(rate_map);
    assert(rate_container.outgoingRateCount(2)==0);
    assert(rate_container.outgoingRateCount(1)==1);
    assert(rate_container.outgoingRateCount(3)==1);
  }

  cout << "Testing: addRate" << endl; 
  {
    double rate = 1.0;
    double rate2 = 2.0;

    KMC_Rate_Container rate_container;
    rate_container.addRate(1,2,&rate);
    rate_container.addRate(3,2,&rate2);

    assert(rate_container.incomingRateCount(2)==2);
    assert(rate_container.incomingRateCount(1)==0);
    assert(rate_container.incomingRateCount(3)==0);
  }

  cout << "Testing: addRates" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;
    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;

    KMC_Rate_Container rate_container;
    rate_container.addRates(rate_map);
  
    int value = rate_container.outgoingRateCount(2);
    assert(value==0);
    value = rate_container.outgoingRateCount(1);
    assert(value==1);
    value = rate_container.outgoingRateCount(3);
    assert(value==1);
  }

  cout << "Testing: getRate" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;

    KMC_Rate_Container rate_container(rate_map);

    const double * rate_ = rate_container.getRate(1,2);
    const double * rate2_ = rate_container.getRate(3,2);

    assert(*rate_==rate);
    assert(*rate2_==rate2);

    bool throw_error = false;
    try{
      rate_container.getRate(2,3);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);

    throw_error = false;
    try{
      rate_container.getRate(0,1);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: getIncomingRates" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;

    KMC_Rate_Container rate_container(rate_map);

    auto incoming_rates = rate_container.getIncomingRates(2);
    assert(incoming_rates.size()==2);

    bool rate_found = false;
    bool rate2_found = false;

    for ( auto rate_it : incoming_rates){
      if( rate_it.first == 1 && rate_it.second.begin()->first == 2) {
        rate_found = true;
        assert(*(incoming_rates[1][2])==rate);
      }
      if( rate_it.first == 3 && rate_it.second.begin()->first == 2) {
        rate2_found = true;
        assert(*(incoming_rates[3][2])==rate2);
      }
    }
    assert(rate_found);
    assert(rate2_found);

    auto incoming_rates2 = rate_container.getIncomingRates(1);
    assert(incoming_rates2.size()==0);
  }

  cout << "Testing: getOutgoingRates" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;

    KMC_Rate_Container rate_container(rate_map);

    auto outgoing_rates = rate_container.getOutgoingRates(2);
    assert(outgoing_rates.size()==0);

    auto outgoing_rates2 = rate_container.getOutgoingRates(1);
    assert(outgoing_rates2.size()==1);
    assert(*outgoing_rates2[1][2]==rate);

    auto outgoing_rates3 = rate_container.getOutgoingRates(3);
    assert(outgoing_rates3.size()==1);
    assert(*outgoing_rates3[3][2]==rate2);
  }

  cout << "Testing: getSourceSiteIds" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;
    double rate3 = 3.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;
    rate_map[4][3] = &rate3;

    KMC_Rate_Container rate_container(rate_map);

    auto sources = rate_container.getSourceSiteIds();

    assert(sources.size()==2);

    bool found_1 = false;
    bool found_4 = false;
    
    for( auto siteId : sources ){
      if(siteId==1) found_1 = true;
      if(siteId==4) found_4 = true;
    }

    assert(found_1);
    assert(found_4);
  }

  cout << "Testing: getSinkSiteIds" << endl;
  {
    double rate = 1.0;
    double rate2 = 2.0;
    double rate3 = 3.0;

    Rate_Map rate_map;
    rate_map[1][2] = &rate; 
    rate_map[3][2] = &rate2;
    rate_map[4][3] = &rate3;

    KMC_Rate_Container rate_container(rate_map);

    auto sinks = rate_container.getSinkSiteIds();
    
    assert(sinks.size()==1);
    assert(sinks.at(0)==2);

  }
  return 0;
}

