#include <catch2/catch.hpp>

#include <cassert>
#include <iostream>

#include "../../libmythical/site_container.hpp"

using namespace std;
using namespace mythical;

TEST_CASE("Testing: Site Container","[unit]"){

  cout << "Testing: Constructor" << endl;
  {
    Site_Container site_container;
  }

  cout << "Testing: addSite" << endl;
  {
    Site_Container site_container;
    
    Site site;
    site.setId(1);

    site_container.addSite(site);

    bool throw_error = false;
    try {
      site_container.addSite(site);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: addSites" << endl;
  {
    Site site;
    site.setId(1);

    vector<Site> sites;
    sites.push_back(site);

    Site_Container site_container;
    site_container.addSites(sites);

    // Attempt to add the same site twice
    sites.push_back(site);
    Site_Container site_container2;
    
    bool throw_error = false;
    try{
      site_container2.addSites(sites);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: getSite" << endl;
  {
    Site site;
    site.setId(1);

    Site_Container site_container;

    site_container.addSite(site);

    auto site2 = site_container.getSite(1);
    assert(site2.getId()==1);

    // Try to grab a site that is not stored in the container
    bool throw_error = false;
    try { 
      site_container.getSite(0);
    }catch(...) {
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: getSites" << endl;
  {
    Site site;
    Site site2;

    site.setId(1);
    site2.setId(2);

    Site_Container site_container;

    site_container.addSite(site);
    site_container.addSite(site2);

    vector<int> siteIds = { 1, 2 };
    auto sites = site_container.getSites(siteIds);

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
      site_container.getSites(siteIds2);
    }catch(...) {
      throw_error = true;
    }
    assert(throw_error);

    // Grab sites without providing arguments should grab all the sites
    sites = site_container.getSites();

    assert(sites.size()==2);

    found1 = false;
    found2 = false;

    for( auto site : sites ){
      if(site.second.getId()==1) found1 = true;
      if(site.second.getId()==2) found2 = true;
    }

    assert(found1);
    assert(found2);
  }

  cout << "Testing: size" << endl;
  {
    Site site;
    Site site2;

    site.setId(1);
    site2.setId(2);

    Site_Container site_container;
    assert(site_container.size()==0);

    site_container.addSite(site);
    site_container.addSite(site2);
    assert(site_container.size()==2);
  }
 
  cout << "Testing: getClusterIdOfSite" << endl;
  {
    Site site;
    site.setId(1);
    site.setClusterId(3);

    Site_Container site_container;
    site_container.addSite(site);

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
    Site site;
    Site site2;

    site.setId(1);
    site2.setId(2);

    site.setClusterId(3);

    Site_Container site_container;
    site_container.addSite(site);
    site_container.addSite(site2);

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
    Site site;
    Site site2;

    site.setId(1);
    site2.setId(2);

    site.setClusterId(3);
    site2.setClusterId(6);

    Site_Container site_container;
    site_container.addSite(site);
    site_container.addSite(site2);

    vector<int> siteIds = { 1 , 2};
    assert(site_container.getSmallestClusterId(siteIds)==3);

  } 

  cout << "Testing: getRateToNeighborOfSite" << endl;
  {
    Site site;
    Site site2;

    site.setId(1);
    site2.setId(2);

    double rate1_2 = 1.0;
    double rate2_1 = 2.0;

    site.addNeighRate(pair<int,double *>(2,&rate1_2));
    site2.addNeighRate(pair<int,double *>(1,&rate2_1));

    Site_Container site_container;
    site_container.addSite(site);
    site_container.addSite(site2);

    assert(site_container.getRateToNeighborOfSite(1,2)==rate1_2);
    assert(site_container.getRateToNeighborOfSite(2,1)==rate2_1);

  }
}
