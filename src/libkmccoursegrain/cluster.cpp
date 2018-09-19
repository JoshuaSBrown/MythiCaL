
#include <cmath>

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
  iterations_ = 3;
  convergenceTolerance_ = 0.01;
  convergence_method_ = converge_by_iterations_per_site;
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

double Cluster::getProbabilityOfOccupyingInternalSite(int siteId){
  if(!sitesInCluster_.count(siteId)){
    throw invalid_argument("the provided site is not in the cluster");
  }
  return probabilityOnSite_[siteId];
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
void Cluster::initializeProbabilityOnSites_(){
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

*/
// int in the map is the neigh 
// int in the vector is the site 
// vector contains rates as shown
//  
// neigh1 -> site1 <- neigh2 
//
//
map<const int,vector<pair<const int, double>>> 
Cluster::getInternalRatesFromNeighborsComingToSite_(){

  map<const int, vector<pair<const int,double>>> internalRates;

  for(auto site : sitesInCluster_){
    for(auto neighId : site.second->getNeighborIds()) {
      if(siteIsInCluster(neighId)){
        pair<const int,double> rateToSite(site.first,site.second->getRateToNeighbor(neighId));
        internalRates[neighId].push_back(rateToSite);
      }
    }
  } 
  return internalRates;
}

double Cluster::getProbabilityOfHoppingToNeighbor(int neighId){
  if(!probabilityHopToNeighbor_.count(neighId)){
    string err = "Cannot get probability of hopping to neighbor, either the"
      " site is not a neighbor or the cluster has not been converged.";
    throw invalid_argument(err);
  }
  return probabilityHopToNeighbor_[neighId];
}

// The firsts int is the id of site within cluster
// second int is the id of the neighbor of site outside of cluster
map<const int, map<const int, double>>
Cluster::getRatesToNeighborsOfCluster_(){

  map<const int, map<const int,double>> externalRates;

  for(auto site : sitesInCluster_){
    for(auto neighId : site.second->getNeighborIds()) {
      if(!siteIsInCluster(neighId)){
        externalRates[site.first][neighId]=site.second->getRateToNeighbor(neighId);
      }
    }
  } 
  return externalRates;
}

void Cluster::setConvergenceTolerance(double tolerance){
  if(tolerance<0.0){
    throw invalid_argument("tolerance must be a positive value");
  }
  convergenceTolerance_ = tolerance;
}

void Cluster::setConvergenceIterations(long iterations){
  if(iterations<1){
    throw invalid_argument("number of iterations must be greater than 0.");
  }
  iterations_ = iterations;
}


void Cluster::iterate_(map<const int, vector<pair<const int,double >>> ratesBetween){

  map<const int, double> temp_probabilityOnSite;

  double total = 0;
  for(auto site : sitesInCluster_){
    for( auto neighsite : site.second->getNeighborIds()){
      if(siteIsInCluster(neighsite)){
        if(temp_probabilityOnSite.count(site.first)){
          temp_probabilityOnSite[site.first] += sitesInCluster_[neighsite]->probHopToNeigh(site.first)*
            probabilityOnSite_[neighsite];
        }else{  
          temp_probabilityOnSite[site.first] = sitesInCluster_[neighsite]->probHopToNeigh(site.first)*
            probabilityOnSite_[neighsite];
        }
        total+=temp_probabilityOnSite[site.first];
      }
    }
  }

  double total2 = 0.0;
  for(auto site : sitesInCluster_){
    probabilityOnSite_[site.first] = (temp_probabilityOnSite[site.first]/total+probabilityOnSite_[site.first])/2.0;
    total2 += probabilityOnSite_[site.first];
  }

  for(auto site: sitesInCluster_){
    probabilityOnSite_[site.first] = probabilityOnSite_[site.first]/total2;
  }

}

void Cluster::calculateProbabilityHopToNeighbors_(){
  
  auto ratesToNeighbors = getRatesToNeighborsOfCluster_();

  cout << "Rates to Neighbors" << endl;


  auto sumRatesOffCluster = 0.0;
  for( auto rateToNeigh : ratesToNeighbors ) {
    for( auto rate : rateToNeigh.second ){
      cout << "from " << rateToNeigh.first << " to " << rate.first << " rate " << rate.second << endl;
      sumRatesOffCluster+=rate.second;
    }
  }

  auto sumDwell = 0.0;
  for( auto site : sitesInCluster_ ){
    sumDwell += site.second->getDwellTime();
  }
  cout << "sumDwell " << sumDwell << endl;

  map<const int, double> probabilityHopToNeighbor;
  double total = 0.0;
  for(auto rateToNeigh : ratesToNeighbors){
    int siteHoppingFrom = rateToNeigh.first;
    cout << "Site hopping from " << siteHoppingFrom << endl;
    for( auto rate : rateToNeigh.second ){
      int siteHoppingTo = rate.first;
      cout << "Site hopping to " << siteHoppingTo << endl;
      if(probabilityHopToNeighbor.count(siteHoppingTo)){
        cout << "Exists" << endl;
        probabilityHopToNeighbor[siteHoppingTo] +=\
          probabilityOnSite_[siteHoppingFrom]*\
          sitesInCluster_[siteHoppingFrom]->getDwellTime()/\
          sumDwell*\
          sitesInCluster_[siteHoppingFrom]->probHopToNeigh(siteHoppingTo)/\
          sumRatesOffCluster;
      }else{

        cout << "Prob on site " << probabilityOnSite_[siteHoppingFrom] << endl;
        cout << "dwell time " << sitesInCluster_[siteHoppingFrom]->getDwellTime() << endl;
        cout << "Probability Hop to any neighbor " << sitesInCluster_[siteHoppingFrom]->probHopToNeigh(siteHoppingTo) << endl;
        cout << "sumRatesOff " << sumRatesOffCluster << endl;
        probabilityHopToNeighbor[siteHoppingTo] =\
          probabilityOnSite_[siteHoppingFrom]*\
          sitesInCluster_[siteHoppingFrom]->getDwellTime()/\
          sumDwell*\
          sitesInCluster_[siteHoppingFrom]->probHopToNeigh(siteHoppingTo)/\
          sumRatesOffCluster;
      }
      total += probabilityHopToNeighbor[siteHoppingTo];
    }
  }

  for(auto neighborProb : probabilityHopToNeighbor ){
    probabilityHopToNeighbor[neighborProb.first] = neighborProb.second/total;
  }
  probabilityHopToNeighbor_ = probabilityHopToNeighbor;
  cout << "Number of elements " << probabilityHopToNeighbor.size() << endl;
  for(auto val : probabilityHopToNeighbor_){
    cout << val.first << " " << val.second << endl;
  }
}

void Cluster::converge(){

  initializeProbabilityOnSites_();
  auto ratesBetweenInternalSites = getInternalRatesFromNeighborsComingToSite_();

  if(convergence_method_==converge_by_iterations_per_cluster){
    for(long i = 0; i < iterations_; i++){
      iterate_(ratesBetweenInternalSites); 
    }
  }else if(convergence_method_==converge_by_iterations_per_site){
    long total_iterations = iterations_*static_cast<long>(sitesInCluster_.size());
    for(long i = 0; i < total_iterations; i++){
      iterate_(ratesBetweenInternalSites); 
    }
  }else{
    double error = convergenceTolerance_*1.1;
    while(error>convergenceTolerance_){
      auto oldSiteProbs = probabilityOnSite_;
      iterate_(ratesBetweenInternalSites); 
      error = 0.0;
      for( auto site : oldSiteProbs ){
        error += pow(oldSiteProbs[site.first]-probabilityOnSite_[site.first],2.0);
      }
      error = pow(error,1.0/2.0);
    }
  }
  calculateProbabilityHopToNeighbors_();
}

}
