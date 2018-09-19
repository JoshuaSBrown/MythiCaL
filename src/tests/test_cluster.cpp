#include <iostream>
#include <cassert>
#include <vector>
#include <memory>
#include <cmath>
//#include <map>
//#include <iterator>

//#include <kmccoursegrain/site.hpp>
#include <kmccoursegrain/cluster.hpp>
#include <kmccoursegrain/site.hpp>

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

  cout << "Testing: convergence1" << endl;
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

    cluster.converge();

    assert(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(1)))==50);
    assert(round(static_cast<int>(100*cluster.getProbabilityOfOccupyingInternalSite(2)))==50);

  }

  cout << "Testing: convergence2" << endl;
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

    cluster.converge();

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

    cluster.converge();

    assert(round(static_cast<int>(100*cluster.getProbabilityOfHoppingToNeighbor(4)))==100);

    double rate6 = 1;
    siteSmart->addNeighRate(pair<int const, double *>(5,&rate6));
    cluster.converge();
    
    // 
    // neigh5 <- site1 -> site2  -> site3 -> neigh4
    //                 <-        <-

    assert(round(static_cast<int>(100*cluster.getProbabilityOfHoppingToNeighbor(4))) ==50);
    assert(round(static_cast<int>(100*cluster.getProbabilityOfHoppingToNeighbor(5))) ==50);

  }

	//Testing the Constructor
/*	{
		cout<<"Testing the constructor"<<endl;
		shared_ptr<Cluster> tmp = make_shared<Cluster>();


		vector<double> nRates = {1.0,2.0};
		vector<int> nId = {2,3};
		map<int, double> neighbors1 = {
			{2,1.0},
			{3,2.0}
		};
		int vFreq1 = 1;
		int sId1 = 1;

		vector<double> nRates2 = {1.0,5.0};
		vector<int> nId2 = {1,3};
		map<int, double> neighbors2 = {
			{1,1.0},
			{3,5.0}
		};
		int vFreq2 = 2;
		int sId2 = 2;

		vector<double> nRates3 = {2.0,5.0};
		vector<int> nId3 = {1,2};
		map<int, double> neighbors3 = {
			{1,2.0},
			{2,5.0}
		};
		int vFreq3 = 3;
		int sId3 = 3;
		//site site1(sId, nRates, nId, vFreq);
		shared_ptr<site> site1(new Site(sId1,vFreq1));
		shared_ptr<site> site2(new Site(sId2,vFreq2));
		shared_ptr<site> site3(new Site(sId3,vFreq3));
		
		map<shared_ptr<site>,double> linkSite1 = {
			{site2,1.0},
			{site3,2.0}
		};

		map<shared_ptr<site>,double> linkSite2 = {
			{site1,1.0},
			{site3,5.0},
			{site2,10.0}
		};

		site1->addNeighbors(linkSite1);
		site2->addNeighbors(linkSite2);
	
		site1->printInfo();
		site2->printInfo();	

		tmp->addSite(site1);
		
		tmp->addSite(site2);
		tmp->printInfo();
		cout<<"Dwell Time: "<<tmp->dwellTime()<<endl;
		cout<<"Prob Hop from site 1 to site 2: "<<site1->probHop(2)<<endl;
		cout<<"Prob Hop from site 1 to site 2 (overloaded):"<<site1->probHop(site2)<<endl;
	
	//Testing Convergence
	
		tmp->convergence(5);
		tmp->printInfo();
		site1->printInfo();
		site2->printInfo();
	}	*/
	return 0;
}
