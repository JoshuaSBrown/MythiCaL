//#include "cluster.h"

#include <kmccoursegrain/site.hpp>

using namespace std;
using namespace kmccoursegrain;

/*
typedef shared_ptr<Site> sitePtr;

Site::Site(int sId, int vFreq, std::map<std::shared_ptr<Site>,double> nSites){
	siteId=sId;
	visitFreq=vFreq;
	clustTag = -1;
	probOnSite = 0.0;
	this->addNeighbors(nSites);
}

Site::Site(int sId, int vFreq){
	siteId=sId;
	visitFreq=vFreq;
	clustTag = -1;
	probOnSite = 0.0; 
}*/

void Site::setRatesToNeighbors(map<int const, double&> neighRates){
  neighRates_ = neighRates;
}
/*
int site::addNeighbors(map<sitePtr, double> addSites){
	int found;

	for(auto it = addSites.cbegin(); it != addSites.cend(); ++it){
		if(it->first == NULL){
			if(Err) cerr<<"ERROR in addNeighbors: adding a null site"<<endl;
			return 0;
		}
		if(it->first->siteId == this->siteId) continue;
		found = 0;
		for(int j = 0; j < (int) neighSites.size(); j++){
			if(neighSites[j]->siteId == it->first->siteId) found = 1;
		}
		
		if(!found) neighSites.push_back(it->first);
		neighs[it->first->siteId] = it->second;

		found = 0;
		for(int j = 0; j < (int) it->first->neighSites.size(); j++){
			if(it->first->neighSites[j]->siteId == this->siteId) found = 1;
		}
		if(!found) it->first->neighSites.push_back(sitePtr(this));
		it->first->neighs[siteId] = it->second;
	}
	return 1;
}
*/

/*
void site::printInfo(){
	cout<<"Site Id: "<<siteId<<endl;
	cout<<"Cluster tag: "<<clustTag<<endl;
	cout<<"Visit Frequency: "<<visitFreq<<endl;
	cout<<"Prob On Site: "<<probOnSite<<endl;
	cout<<"Neighbors:Rates"<<endl;
	for(auto it = neighs.cbegin(); it != neighs.cend(); ++it){
		cout<<"\t"<<it->first<<":"<<it->second<<endl;
	}
	return;
}

double site::getProbOnSite(){
	return probOnSite;
}
	*/ 

/*
//Overload the probHop function
double site::probHop(int receiving){
	bool errFlag = true;
	double totalRates = 0.0;

	for(auto it = neighs.cbegin(); it != neighs.cend(); ++it){
	       if(it->first == receiving) errFlag = false;
	       totalRates += it->second;
	}

	if(errFlag){
		if(Err) cerr<<"ERROR in probHop: receiving site not a neighbor of shipping site"<<endl;
		return 0;
	}
	return (neighs[receiving]/totalRates);
}

double site::probHop(sitePtr receiving){ 
	bool errFlag = true;
	double totalRates = 0.0;

	for(auto it = neighs.cbegin(); it != neighs.cend(); ++it){
	       if(it->first == receiving->siteId) errFlag = false;
	       totalRates += it->second;
	}

	if(errFlag){
		if(Err) cerr<<"ERROR in probHop: receiving site not a neighbor of shipping site"<<endl;
		return 0;
	}
	return (neighs[receiving->siteId]/totalRates);
}
*/


