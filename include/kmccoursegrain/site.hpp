#ifndef KMCCOURSEGRAIN_CLUSTER_H_
#define KMCCOURSEGRAIN_CLUSTER_H_

#include <vector>
#include <memory>

namespace kmccoursegrain{

struct Site{
	int siteId;
	std::vector<double> neighRates;
	std::vector<int> neighIds; 
	int visitFreq;
	int clustTag;
	double ProbOnSite;

	site(){};//default constructor

	//constructor
	site(int sId, std::vector<double> nRates, std::vector<int> nId, int vFreq){
		//error handling?
		siteId=sId;
		visitFreq=vFreq;
		clustTag = -1;
		ProbOnSite = 0.0;
		//optional way of handeling it
		/*
		for(int i; i<nId.size(); i++){
			neighRates.push_back(nRates[i])
			neighIds.push_back(nId[i]);
			}*/
		neighRates.swap(nRates);
		neighIds.swap(nId);
	}
};

}
#endif // KMCCOURSEGRAIN_CLUSTER_H_
