#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "cluster.h"

using namespace std;

#ifdef _E_
#define Err 1
#else
#define Err 0
#endif

int main(void){
	
	if(Err==1) cout<<"Error Reporting On"<<endl;
	//Testing the Threshold Setter
	{
		cout<<"Testing Thresh"<<endl;
		assert(setThresh(10));
		assert(getThresh()==10);
		assert(setThresh(-20)==0);
	}
	//Testing the Constructor
	{
		cout<<"Testing the constructor"<<endl;
		shared_ptr<cluster> tmp = make_shared<cluster>();
		vector<double> nRates = {1.0,2.0};
		vector<int> nId = {2,3};
		int vFreq = 1;
		int sId = 1;

		vector<double> nRates2 = {1.0,3.0};
		vector<int> nId2 = {1,3};
		int vFreq2 = 2;
		int sId2 = 2;
		//site site1(sId, nRates, nId, vFreq);
		shared_ptr<site> site1(new site(sId,nRates,nId,vFreq));
		shared_ptr<site> site2(new site(sId2,nRates2,nId2,vFreq2));

		tmp->addSite(site1);
		tmp->addSite(site2);
		tmp->printClusterInfo();
		cout<<"Dwell Time: "<<tmp->dwellTime()<<endl;
		cout<<"Prob Hop from site 1 to site 2: "<<tmp->probHop(site1,site2)<<endl;
	}
	return 0;
}
