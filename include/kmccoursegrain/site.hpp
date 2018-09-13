#ifndef KMCCOURSEGRAIN_SITE_H_
#define KMCCOURSEGRAIN_SITE_H_

#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <math.h>

//All functions return 1 if success, 0 or nan if failure as well as cerr printing
struct Site{
//	friend class cluster;
	public:	
		Site() : probabilitySet_(false), visitFreq_(0) {};
    void setRatesToNeighbors(std::map<int const, double&> neighRates);
/*
		//Two ways to create a site object
		//Creates site with list of sites and rates as neighbors
		site(int sId, int vFreq, std::map<std::shared_ptr<site>,double> nSites);
		//Just creates site
		site(int sId, int vFreq);
		//Two ways to call probHop hop is from shipping to receivingSite
		double probHop(std::shared_ptr<site> receivingSite);
		double probHop(int receivingSiteId);
		int addNeighbors(std::map<std::shared_ptr<site>, double>);
		void printInfo();
		double getProbOnSite();
		int mergeSites();//FIXME
*/
	private:
//		int siteId;
//		std::vector<std::shared_ptr<site>> neighSites;
//		std::map<int, double> neighs;
    std::map<int const, double*> neighRates_;
		int visitFreq_;
    double courseGrainedSiteProb_;
    bool probabilitySet_;
//		int clustTag;
//		double probOnSite;
};

#endif
