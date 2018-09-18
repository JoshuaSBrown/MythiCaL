//#include "cluster.h"

#include <kmccoursegrain/site.hpp>

using namespace std;

namespace kmccoursegrain {
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

void Site::setRatesToNeighbors(map<int const, double*> neighRates){
  neighRates_ = neighRates;
}

void Site::addNeighRate(pair<int const, double*> neighRate){
  if(neighRates_.count(neighRate.first)){
    throw invalid_argument("That neighbor has already been added.");
  }
  neighRates_[neighRate.first] = neighRate.second;
}

void Site::resetNeighRate(pair<int const, double*> neighRate){
  neighRates_[neighRate.first] = neighRate.second;
}

vector<int> Site::getNeighborIds(){
  vector<int> neighborIds;
  for(auto neighId : neighRates_) neighborIds.push_back(neighId.first);
  return neighborIds;
}

vector<double> Site::getRatesToNeighbors(){
  vector<double> rates;
  for(auto rate : neighRates_) rates.push_back(*(rate.second));
  return rates;
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

std::ostream& operator<<(std::ostream& os, const kmccoursegrain::Site& site){
  os << "Site Id: "<<site.getId() << std::endl;
  os << "Cluster Id: "<<site.clusterId_ << std::endl;
  os << "Visit Frequency: "<<site.visitFreq_ << std::endl;
  os << "Neighbors:Rates"<<std::endl;
  for(auto rate_ptr : site.neighRates_ ){
    os << "\t"<<rate_ptr.first<<":"<< *(rate_ptr.second) << std::endl;
  }
  return os;
}

//Overload the probHop function
double Site::probHopToNeigh(const int neighSiteId){

  if(neighRates_.count(neighSiteId)==0){
    string err = "Error site " +to_string(neighSiteId)+" is not a nieghbor of "
      "" + to_string(getId());
    throw invalid_argument(err);
  }

  double totalRates = 0.0;
  for(auto ptr_rate : neighRates_) totalRates += *(ptr_rate.second);

	return (*(neighRates_[neighSiteId])/totalRates);
}

/*
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
}
