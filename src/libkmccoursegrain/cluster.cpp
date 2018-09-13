
#include <kmccoursegrain/cluster.hpp>

using namespace std;
using namespace kmccoursegrain;

#ifdef _E_
#define Err 1
#else
#define Err 0
#endif


//typedef shared_ptr<site> sitePtr;

static int clusterIdCounter = 0;
static int thresh = 0;

int setThresh(int n){
    if(n>0){
        thresh=n;
        if(Err) cout<<"Thresh hold is: "<<thresh<<endl;
        return 1;
    }else{
        if(Err) cerr<<"ERROR in setThresh: threshold is negative"<<endl;
        return 0;
    }
}

int getThresh(){
	return thresh;
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
cluster::cluster(){
	clustTag=clusterIdCounter;
	clusterIdCounter++;
}

cluster::~cluster(){
}

int cluster::addSite(sitePtr newSite){
	if(newSite == NULL){
		if(Err) cerr<<"ERROR in addSite: adding a null site"<<endl;
		return 0;
	}

	newSite->clustTag = clustTag;

	sitesInCluster.push_back(newSite);
	
	int addFlag;
	for(auto it = newSite->neighs.cbegin(); it != newSite->neighs.cend(); ++it){
		addFlag = 1;
		for(auto jt = neighs.cbegin(); jt != neighs.cend(); ++jt){
			if(jt->first == newSite->siteId){
				neighs.erase(jt);
			}
			if(it == jt){
				addFlag = 0;
			}
		}
		if(addFlag){
			neighs.insert(it, pair<int, double>(it->first,it->second));
		}
	}
	return 1;
}

int cluster::getClusterId(){
	return clustTag;
}


vector<sitePtr> cluster::getSitesInCluster(){
	return sitesInCluster;
}


void cluster::printInfo(){
	cout<<"Cluster ID: "<<clustTag<<endl;
	cout<<"Visit Frequency to Cluster: "<< visitFreq<<endl;
	cout<<"Sites in Cluster: "<<sitesInCluster.size()<<endl;
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
}

double cluster::dwellTime(){
	sitePtr tmp;
	double sum;
	for(int i =0; i < (int)sitesInCluster.size(); i++){
		tmp=sitesInCluster[i];
		for(auto it = tmp->neighs.cbegin(); it != tmp->neighs.cend(); ++it){
			sum += it->second;
		}
	}
	return (1/sum);
}

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

void cluster::initProbOnSite(){
	for(int i =0; i < (int) sitesInCluster.size(); i++){
		sitesInCluster[i]->probOnSite = 1.0/((double) sitesInCluster.size());
	}
	if(Err) cout<<"Intializing Probs On site"<<endl;
	return;
}

int cluster::convergence(long iterations){
	vector<double> interProbs(sitesInCluster.size());
	double total = 0;
	double norm = 0;
	initProbOnSite();
	
	if(iterations <= 0){
		if(Err) cerr<<"ERROR in convergence: iterations is not valid"<<endl;
		return 0;
	}
	for(long i = 0; i < iterations; i++){
		fill(interProbs.begin(),interProbs.end(), 0);
		for(int j = 0; j < (int)sitesInCluster.size(); j++){
			for(int k = 0; k < (int)sitesInCluster[j]->neighSites.size(); k++){
				if(sitesInCluster[j]->neighSites[k]->clustTag == clustTag){
					interProbs[j] += sitesInCluster[j]->neighSites[k]->probOnSite * 
						sitesInCluster[j]->neighSites[k]->probHop(sitesInCluster[j]);
				}
			}
			total += interProbs[j];	
		}
		for(int j = 0; j < (int) sitesInCluster.size(); j++){
			sitesInCluster[j]->probOnSite = ((sitesInCluster[j]->probOnSite + 
				(interProbs[j]/total)) / 2);
			norm += sitesInCluster[j]->probOnSite;
		}
		for(int j = 0; j < (int) sitesInCluster.size(); j++){
			sitesInCluster[j]->probOnSite = (sitesInCluster[j]->probOnSite / norm);
		}
		total = 0;
		norm = 0;
	}
	return 1;		
}*/

