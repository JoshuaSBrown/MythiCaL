#include <iostream>
#include <cassert>
#include <vector>
#include <memory>
#include <cmath>

#include "../libkmccoursegrain/cluster.hpp"
#include "../libkmccoursegrain/site.hpp"

using namespace std;
using namespace kmccoursegrain;

int main(void){

//  Site site;

  cout << "Testing: Cluster constructor" << endl;
  {
    Cluster cl;
  }

  cout << "Testing: Cluster identity setter" << endl;
  {
    Cluster cl;
    cl.setId(0);
  }

  cout << "Testing: Cluster identity getter" << endl;
  {
    Cluster cl;
    cl.setId(0);
    assert(cl.getId()==0);

    bool fail = false;
    Identity cl2;
    try {
      cl2.getId();
    }catch(...){
      fail = true;
    }
    assert(fail);
  }  

	//Testing the Threshold Setter
  cout << "Testing: Threshold setter and getter" << endl;
	{
		setThreshold(10);
		assert(getThreshold()==10);
	}

  cout << "Testing: addSite" << endl;
  {
    Site site;
    site.setId(1);
    double rate = 1.0;
    site.addNeighRate(pair<int const, double *>(2,&rate));

    Cluster cluster;
    cluster.addSite(make_shared<Site>(site));
    assert(cluster.getNumberOfSitesInCluster()==1);
  }


  cout << "Testing: siteIsInCluster" << endl;
  {
    Site site;
    site.setId(1);
    
    Cluster cluster;
    assert(!cluster.siteIsInCluster(1));
    cluster.addSite(make_shared<Site>(site));
    assert(cluster.siteIsInCluster(1));
  }

  cout << "Testing: getProbabilityOfOccupyingInternalSite 1" << endl;
  {

    // Simple convergence test 
    // 
    // site1 -> site2
    //       <-
    //
    // Same rate should lead to 50 % probability on either site

    Site site;
    site.setId(1);
    double rate = 1;
    site.addNeighRate(pair<int const, double *>(2,&rate));
    
    Site site2;
    site2.setId(2);
    double rate2 = 1;
    site2.addNeighRate(pair<int const, double *>(1,&rate2));
  
    Cluster cluster;
    cluster.addSite(make_shared<Site>(site));
    cluster.addSite(make_shared<Site>(site2));

    assert(static_cast<int>(round(100*cluster.getProbabilityOfOccupyingInternalSite(1)))==50);
    assert(static_cast<int>(round(100*cluster.getProbabilityOfOccupyingInternalSite(2)))==50);

  }

  cout << "Testing: getProbabilityOfOccupyingInternalSite 2" << endl;
  {

    // Simple convergence test 
    // 
    // site1 -> site2  -> site3
    //       <-        <-
    //
    // Same rate should lead to 25 % probability on end sites
    // and 50 % probability on middle site

    Site site;
    site.setId(1);
    double rate = 1;
    site.addNeighRate(pair<int const, double *>(2,&rate));
    
    Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int const, double *>(1,&rate2));
    site2.addNeighRate(pair<int const, double *>(3,&rate3));
  
    Site site3;
    site3.setId(3);
    double rate4 = 1;
    site3.addNeighRate(pair<int const, double *>(2,&rate4));

    Cluster cluster;
    cluster.setConvergenceIterations(6);

    cluster.addSite(make_shared<Site>(site));
    cluster.addSite(make_shared<Site>(site2));
    cluster.addSite(make_shared<Site>(site3));

    assert(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(1)))==25);
    assert(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(2)))==50);
    assert(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(3)))==25);

  }

  cout << "Testing: get probability of hopping to a neighbor" << endl;
  {
    // Simple convergence test 
    // 
    // site1 -> site2  -> site3 -> neigh4
    //       <-        <-
    //
    // Same rate should lead to 25 % probability on end sites
    // and 50 % probability on middle site

    Site site;
    site.setId(1);
    double rate = 1;
    site.addNeighRate(pair<int const, double *>(2,&rate));
    auto siteSmart = make_shared<Site>(site);   
 
    Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int const, double *>(1,&rate2));
    site2.addNeighRate(pair<int const, double *>(3,&rate3));
  
    Site site3;
    site3.setId(3);
    double rate4 = 1;
    site3.addNeighRate(pair<int const, double *>(2,&rate4));
    double rate5 = 1;
    site3.addNeighRate(pair<int const, double *>(4,&rate5));

    Cluster cluster;
    cluster.setConvergenceIterations(6);

    cluster.addSite(siteSmart);
    cluster.addSite(make_shared<Site>(site2));
    cluster.addSite(make_shared<Site>(site3));

    assert(round(static_cast<int>(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(4)))==100);

    double rate6 = 1;
    siteSmart->addNeighRate(pair<int const, double *>(5,&rate6));
   
    cluster.updateProbabilitiesAndTimeConstant(); 
    // 
    // neigh5 <- site1 -> site2  -> site3 -> neigh4
    //                 <-        <-
    assert(static_cast<int>(round(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(4))) ==50);
    assert(static_cast<int>(round(100*cluster.getProbabilityOfHoppingToNeighborOfCluster(5))) ==50);

  }

  cout << "Testing: pickNewSiteId" << endl;
  {
    Site site;
    site.setId(1);
    double rate = 1;
    double rate6 = 1;
    site.addNeighRate(pair<int const, double *>(2,&rate));
    site.addNeighRate(pair<int const, double *>(5,&rate6));
 
    Site site2;
    site2.setId(2);
    double rate2 = 1;
    double rate3 = 1;
    site2.addNeighRate(pair<int const, double *>(1,&rate2));
    site2.addNeighRate(pair<int const, double *>(3,&rate3));
  
    Site site3;
    site3.setId(3);
    double rate4 = 1;
    site3.addNeighRate(pair<int const, double *>(2,&rate4));
    double rate5 = 1;
    site3.addNeighRate(pair<int const, double *>(4,&rate5));

    Cluster cluster;
    cluster.setConvergenceIterations(6);

    cluster.addSite(make_shared<Site>(site));
    cluster.addSite(make_shared<Site>(site2));
    cluster.addSite(make_shared<Site>(site3));

    // Setting the seed will ensure that the results are reproducable

    cluster.setRandomSeed(1);
    
    assert(cluster.pickNewSiteId()==5); 
    assert(cluster.pickNewSiteId()==3); 
    assert(cluster.pickNewSiteId()==2); 
    assert(cluster.pickNewSiteId()==2); 
    assert(cluster.pickNewSiteId()==3); 
    assert(cluster.pickNewSiteId()==2); 
    assert(cluster.pickNewSiteId()==1); 
    assert(cluster.pickNewSiteId()==3); 
  }

	return 0;
}
