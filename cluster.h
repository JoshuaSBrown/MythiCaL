#ifndef CLUSTER_H_
#define CLUSTER_H_
#include <vector>
#include <memory>
struct site{
	int siteId;
	std::vector<double> neighRates;
	std::vector<int> neighIds; 
	int visitFreq;
	site(){};//default constructor
	site(int sId, std::vector<double> nRates, std::vector<int> nId, int vFreq){
		siteId=sId;
		visitFreq=vFreq;
		neighRates.swap(nRates);
		neighIds.swap(nId);
	}
};
class cluster {
	public:
		cluster();
		int addSite(std::shared_ptr<site>); //addSites to cluster
		int getClusterId();
		int printClusterInfo();
		double dwellTime();
	private:
		int clusterId; //should also be included in site struct
		int visitFreqCluster;
};
#endif
