
#include <kmccoursegrain/cluster.hpp>
#include <kmccoursegrain/site.hpp>
#include <kmccoursegrain/log.hpp>

using namespace std;

namespace kmccoursegrain {

typedef shared_ptr<Site> SitePtr;

static int clusterIdCounter = 0;
static int threshold = 0;

void setThreshold(const int n){
    if(n>0){
        threshold=n;
        LOG("Threshold is "+to_string(threshold),1);
    }else{
        throw runtime_error("ERROR in setThresh, threshold is negative");
    }
}

int getThreshold(){
	return threshold;
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


Cluster::Cluster(){
	setId(clusterIdCounter);
	clusterIdCounter++;
}


void Cluster::addSite(SitePtr newSite){
	if(newSite == NULL){
		throw invalid_argument("ERROR in addSite: adding a null site");
	}
  if(sitesInCluster_.count(newSite->getId())){
    throw invalid_argument("Site has already been added to the cluster");
  }
	newSite->setClusterId(this->getId());
	sitesInCluster_[newSite->getId()]=newSite;
}

vector<SitePtr> Cluster::getSitesInCluster(){
  vector<SitePtr> sites;
  for( auto site : sitesInCluster_ ) sites.push_back(site.second);
	return sites;
}

std::ostream& operator<<(std::ostream& os, const kmccoursegrain::Cluster& cluster){
  os << "Cluster Id: "<<cluster.getId() << std::endl;
  os << "Cluster visitFreq: "<<cluster.visitFreqCluster_ << std::endl;
	os << "Number of sites in Cluster: " << cluster.sitesInCluster_.size()<<endl;
  os << "Sites in cluster: " << endl;
  for( auto site : cluster.sitesInCluster_ ){
    os << *(site.second) << endl;
  }
  return os;
}
/*
void cluster::printInfo(){
	cout<<"Cluster ID: "<<clustTag<<endl;
	cout<<"Visit Frequency to Cluster: "<< visitFreq<<endl;
	for(int i = 0; i < (int)sitesInCluster.size(); i++){
		cout<<"Site Id: "<<sitesInCluster[i]->siteId<<" Visit Frequency: "<<sitesInCluster[i]->visitFreq<<endl;
	}
	std::map<int, double> neighbors;
	cout<<"Rates: "<<endl;
	for(int i = 0; i < (int)sitesInCluster.size(); i++){
		cout<<"\tSite: "<<sitesInCluster[i]->siteId<<endl;
		neighbors = sitesInCluster[i]->neighs;
		cout<<"\t\tNeighbor:Rate"<<endl;
		for(auto it = neighbors.cbegin(); it != neighbors.cend(); ++it){
			cout<<"\t\t"<<it->first<<":"<<it->second<<endl;
		}
	}
	return;
}*/

double Cluster::dwellTime(){

  throw runtime_error("Dwelltime needs to be fixed current implementation is not correct.");
	SitePtr site_ptr;
	double sum = 0.0;
	for(auto site : sitesInCluster_){ 
		for(auto rate : site.second->getRatesToNeighbors()){
			sum += rate;
		}
	}
	return (1/sum);
}
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


double cluster::probHopOff(sitePtr target, long iterations){
	if(iterations <= 0){
		if(Err) cerr<<"ERROR in probHopOff: iterations is not valid"<<endl;
		return nan("");
	}		
	convergence(iterations);
	double interSum = 0.0;
	double sum = 0.0;
	double targetSum = 0.0;
	double tmpInter = 0.0;
	
	if(target == NULL){
		if(Err) cerr<<"ERROR in probHopOff: targeting a null site"<<endl;
		return 0;
	}

	int found = 1;
	for(int i = 0; i < (int) sitesInCluster.size(); i++){
		for(int j = 0; j < (int) sitesInCluster[i]->neighSites.size(); j++){
			if(sitesInCluster[i]->neighSites[j]->clustTag != clustTag){
				tmpInter = sitesInCluster[i]->probHop(sitesInCluster[i]->neighSites[j]);
				interSum += tmpInter;
				if(sitesInCluster[i]->neighSites[j]->siteId == target->siteId){
					found = 0;
					targetSum = tmpInter * sitesInCluster[i]->probOnSite;
				}
			}
		}
		interSum *= sitesInCluster[i]->probOnSite;
		sum += interSum;
		interSum = 0.0;
	}
	if(!found){
		if(Err) cerr<<"ERROR in probHopOff: target not connected to cluster"<<endl;
		return nan("");
	}
	return(targetSum/sum); 
}

//Overload probHopOff
double cluster::probHopOff(int target, long iterations){
	if(iterations <= 0){
		if(Err) cerr<<"ERROR in probHopOff: iterations is not valid"<<endl;
		return nan("");
	}
	convergence(iterations);
	double interSum = 0.0;
	double sum = 0.0;
	double targetSum = 0.0;
	double tmpInter = 0.0;
	
	int found = 1;
	for(int i = 0; i < (int) sitesInCluster.size(); i++){
		for(int j = 0; j < (int) sitesInCluster[i]->neighSites.size(); j++){
			if(sitesInCluster[i]->neighSites[j]->clustTag != clustTag){
				tmpInter = sitesInCluster[i]->probHop(sitesInCluster[i]->neighSites[j]);
				interSum += tmpInter;
				if(sitesInCluster[i]->neighSites[j]->siteId == target){
					targetSum = tmpInter * sitesInCluster[i]->probOnSite;
				}
			}
		}
		interSum *= sitesInCluster[i]->probOnSite;
		sum += interSum;
		interSum = 0.0;
	}
	if(!found){
		if(Err) cerr<<"ERROR in probHopOff: target not connected to cluster"<<endl;
		return nan("");
	}
	return(targetSum/sum); 
}
*/
void Cluster::initializeProbabilityOnSite_(){
  for(auto site : sitesInCluster_ ){
		probabilityOnSite_[site.first] = 1.0/(static_cast<double>(sitesInCluster_.size()));
	}
	return;
}
/*

// int in the map is the site 
// int in the vector is the neighbor 
// vector contains rates as shown
//  
// neigh1 <- site1 -> neigh2 
//
//
map<const int,vector<pair<const int, double>>> 
Cluster::getInternalRatesFromSiteGoingToNeighbor(){

  map<const int, vector<pair<const int,double>>> internalRates;

  set<int> neighIds;
  for(auto site : sitesInCluster_) neighIds.insert(site.first);

  for(auto site : sitesInCluster_){
    for(auto neighId : site->getNeighIds()) {
      if(siteIsInCluster(neighId)){
        pair<const int,double> rateToNeighbor(neighId,site->getRateToNeighbor(neighId));
        internalRates[site.first].push_back(rateToNeighbor);
      }
    }
  } 
  return internalRates;
}

// int in the map is the neigh 
// int int the vector is the site 
// vector contains rates as shown
//  
// neigh1 -> site1 <- neigh2 
//
//
map<const int,vector<pair<const int, double>>> 
Cluster::getInternalRatesFromNeighborsComingToSite(){

  map<const int, vector<pair<const int,double>>> internalRates;

  set<int> neighIds;
  for(auto site : sitesInCluster_) neighIds.insert(site.first);

  for(auto site : sitesInCluster_){
    for(auto neighId : site->getNeighIds()) {
      if(siteIsInCluster(neighId)){
        pair<const int,double> rateToSite(site.first,site->getRateToNeighbor(neighId));
        internalRates[neighId].push_back(rateToSite);
      }
    }
  } 
  return internalRates;
}


int Cluster::convergence(long iterations){
	double total = 0;
	double norm = 0;
	initProbOnSite();
	
	if(iterations <= 0){
		if(Err) cerr<<"ERROR in convergence: iterations is not valid"<<endl;
		return 0;
	}

	map<int,double> internalSiteProbability;
  double prob = 1.0/static_cast<double>(sitesInCluster_.size());
  for(auto site : sitesInCluster_) internalSiteProbability[site.first] = prob;

  auto ratesToInternalSites = FromNeighborsComingToSite();

	for(long i = 0; i < iterations; i++){

//		for(int j = 0; j < (int)sitesInCluster.size(); j++){
      for(auto site : sitesInCluster_){
//			for(int k = 0; k < (int)sitesInCluster[j]->neighSites.size(); k++){
//        for( auto neighsite : site.getIdsOfNeighSites()){
//				if(sitesInCluster[j]->neighSites[k]->clustTag == clustTag){
  
        internalSiteProbability[site->getId()] +=         
//					internalSiteProbability[j] += sitesInCluster[j]->neighSites[k]->probOnSite * 
						sitesInCluster[j]->neighSites[k]->probHop(sitesInCluster[j]);
				}
			}
			total += internalSiteProbability[j];	
		}
		for(int j = 0; j < (int) sitesInCluster.size(); j++){
			sitesInCluster[j]->probOnSite = ((sitesInCluster[j]->probOnSite + 
				(internalSiteProbability[j]/total)) / 2);
			norm += sitesInCluster[j]->probOnSite;
		}
		for(int j = 0; j < (int) sitesInCluster.size(); j++){
			sitesInCluster[j]->probOnSite = (sitesInCluster[j]->probOnSite / norm);
		}
		total = 0;
		norm = 0;
	}
	return 1;		
}
*/
}
