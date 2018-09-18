#include <iostream>
#include <cassert>
#include <vector>
#include <memory>
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
		cout<<"Testing Thresh"<<endl;
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
