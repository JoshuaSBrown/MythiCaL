
#include <cassert>
#include <iostream>

#include "../../libkmccoarsegrain/kmc_site_container.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

  cout << "Testing: Constructor" << endl;
  {
    KMC_Site_Container site_container;
  }

  cout << "Testing: addKMC_Site" << endl;
  {
    KMC_Site_Container site_container;
    
    KMC_Site site;
    site.setId(1);

    site_container.addKMC_Site(site);

    bool throw_error = false;
    try {
      site_container.addKMC_Site(site);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: addKMC_Sites" << endl;
  {
    KMC_Site site;
    site.setId(1);

    vector<KMC_Site> sites;
    sites.push_back(site);

    KMC_Site_Container site_container;
    site_container.addKMC_Sites(sites);

    // Attempt to add the same site twice
    sites.push_back(site);
    KMC_Site_Container site_container2;
    
    bool throw_error = false;
    try{
      site_container2.addKMC_Sites(sites);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: getKMC_Site" << endl;
  {
    KMC_Site site;
    site.setId(1);

    KMC_Site_Container site_container;

    site_container.addKMC_Site(site);

    auto site2 = site_container.getKMC_Site(1);
    assert(site2.getId()==1);

    // Try to grab a site that is not stored in the container
    bool throw_error = false;
    try { 
      site_container.getKMC_Site(0);
    }catch(...) {
      throw_error = true;
    }
    assert(throw_error);
  }
/*
  cout << "Testing: getKMC_Sites" << endl;
  {
    KMC_Site site;
    KMC_Site site2;

    site.setId(1);
    site2.setId(2);

    KMC_Site_Container site_container;

    site_container.addKMC_Site(site);
    site_container.addKMC_Site(site2);

    vector<int> siteIds = { 1, 2 };
    auto sites = site_container.getKMC_Sites(siteIds);

    assert(sites.size()==2);

    bool found1 = false;
    bool found2 = false;

    for( auto site : sites ){
      if(site.second.getId()==1) found1 = true;
      if(site.second.getId()==2) found2 = true;
    }

    assert(found1);
    assert(found2);

    // Try to find at least one site that is not in the container
    vector<int> siteIds2 = { 1, 3 };

    bool throw_error = false;
    try{
      site_container.getKMC_Sites(siteIds2);
    }catch(...) {
      throw_error = true;
    }
    assert(throw_error);

    // Grab sites without providing arguments should grab all the sites
    sites = site_container.getKMC_Sites();

    assert(sites.size()==2);

    found1 = false;
    found2 = false;

    for( auto site : sites ){
      if(site.second.getId()==1) found1 = true;
      if(site.second.getId()==2) found2 = true;
    }

    assert(found1);
    assert(found2);
  }*/

  cout << "Testing: size" << endl;
  {
    KMC_Site site;
    KMC_Site site2;

    site.setId(1);
    site2.setId(2);

    KMC_Site_Container site_container;
    assert(site_container.size()==0);

    site_container.addKMC_Site(site);
    site_container.addKMC_Site(site2);
    assert(site_container.size()==2);
  }
 
  cout << "Testing: getClusterIdOfSite" << endl;
  {
    KMC_Site site;
    site.setId(1);
    site.setClusterId(3);

    KMC_Site_Container site_container;
    site_container.addKMC_Site(site);

    assert(site_container.getClusterIdOfSite(1)==3);

    bool throw_error = false;
    try{
      site_container.getClusterIdOfSite(0);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: partOfCluster" << endl;
  {
    KMC_Site site;
    KMC_Site site2;

    site.setId(1);
    site2.setId(2);

    site.setClusterId(3);

    KMC_Site_Container site_container;
    site_container.addKMC_Site(site);
    site_container.addKMC_Site(site2);

    bool value = site_container.partOfCluster(1);
    assert(value);
    value =site_container.partOfCluster(2);
    assert(value==false);

    bool throw_error = false;
    try{
      site_container.partOfCluster(0);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);

  }

  cout << "Testing: getSmallestClusterId" << endl;
  {
    KMC_Site site;
    KMC_Site site2;

    site.setId(1);
    site2.setId(2);

    site.setClusterId(3);
    site2.setClusterId(6);

    KMC_Site_Container site_container;
    site_container.addKMC_Site(site);
    site_container.addKMC_Site(site2);

    vector<int> siteIds = { 1 , 2};
    assert(site_container.getSmallestClusterId(siteIds)==3);

  } 

  cout << "Testing: getRateToNeighborOfSite" << endl;
  {
    KMC_Site site;
    KMC_Site site2;

    site.setId(1);
    site2.setId(2);

    double rate1_2 = 1.0;
    double rate2_1 = 2.0;

    unordered_map<int,double> rates;
    rates[2] = rate1_2;
    site.setRatesToNeighbors(&rates);
    unordered_map<int,double> rates2;
    rates2[1] = rate2_1;
    site2.setRatesToNeighbors(&rates2);

    KMC_Site_Container site_container;
    site_container.addKMC_Site(site);
    site_container.addKMC_Site(site2);

    assert(site_container.getRateToNeighborOfSite(1,2)==rate1_2);
    assert(site_container.getRateToNeighborOfSite(2,1)==rate2_1);

  }
  return 0;
}
