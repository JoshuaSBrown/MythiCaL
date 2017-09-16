#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <vector>
#include <memory>

using namespace std;

//functions needed:
//dwell time
//prob hop off to 

struct site{
	int siteId;
	vector<double> neighRates;
	vector<int> neighIds; 
	int visitFreq;

	site(){};//default constructor

	//constructor
	site(int sId, vector<double> nRates, vector<int> nId, int vFreq){
		//error handling?
		siteId=sId;
		visitFreq=vFreq;
		/*
		for(int i; i<nId.size(); i++){
			neighRates.push_back(nRates[i])
			neighIds.push_back(nId[i]);
		} */
		neighRates.swap(nRates);
		neighIds.swap(nId);
	}
};


class cluster {
	public:
		cluster(int            siteId1
			vector<int>    neighIds1
			vector<double> neighRates1
			int            visitFreq1
			int            siteId2
			vector<int>    neighIds2
			vecotr<double> neighRates2
			int            visitFreq2); //constructor  passed certain values? site ids the rates
		~cluster(); //deconstructor
		//add site to cluster function
		int getClusterId();
		//isite ** getSitesInCluster();
		int printClusterInfo();
		// test function and simulation fuction
		double dwellTime();
		//returns 1/sum of rates to neighbors in the cluster

//		int testCluster(void);
//runs all test cases of every function, prints out correct output and what the function prints like such f: (correct) real
//any errors will be reported in the error stream with the funciton call

	private:
		int clusterId; //should also be included in site struct
		int visitFreqCluster;
		vector<shared_ptr<site>> sitesInCluster;
		//add list of neighbors, number of neighbors
		//include an add to cluster funtion to add neighbors of site

//calculate site ratio given hop rates off site, need list of neighbors for that site, hop rates to the neigh form site, list of sites Ids in cluster
//pass maybe cluster struct and return array of site ratios of the sites hop off/hop on see matlab file
//
//fn 2 determine sites in cluster and return list of id sites in cluster,
//
//struct cluster of 2 lsits ids in cluster and neighbors, int id of cluster, int # of sites in cluster, int # of neighbors
//
//create a repository and error handeling
//
//dwell time 1/sum of rates; list of neighbor rates; return double
//
//calculate pvals for site rate(#)/sum of rates site #; list of neighbor hop rates, out put an array
//do for all sites in and out the cluster
//
//calculate p hop off given lisst of neigh rates for a site, output an array
//
//
//Look here 
//dwell time = 1/(sum of rates to neighbors)
//site prob hop = rate to specific site / (sum of rates to all neighbors)
//prob to hop off cluster = rate to hop to neighbor from cluster / (sum of rates to all neighbors to cluster)
//NOTE: prob hop off cluster to site in cluster is 0
//prob on site INITIALLY is 1/number of sites in cluster
//run convergence to fix prob on site in cluster
//t escape = 1/(rates to escape from site out of cluster)
//prob off cluster on neighbor = prob on site in cluster * (dwell time of site in cluster/total dwell time of all sites in cluster) * (rate off to neighbor from site/ total rate off cluster to neighbors)
//
//Return prob on site in cluster prob off cluster on neighbor
//
//convergence
//	1. prob on site is 1/(number of sites in cluster) (only intially)
//	2. site 1 eg (prob hope to site 1)*(prob on site 2)+(prob hppe to site 1)(prob on site 5)
//	3.Total (new site prob)     total=site1prob+site2prob...
//	4.probsite1=site1prob/total + probsite1 ....
//	5.normalise probonsite1=probonsite1/sum(probs on site #)
//
//calculate total prob off cluster to a given neighbor n
//	1. sum hop rates off the cluster
//	2. calc prob hop to neigh off cluster ex for neigh1 prob neigh=dwellofsite(1)/sum(dwell + probonsite1*(rate from site1 to neigh1/tot    tot =sum of hop rates off
//
//is site(1) above refer to as the cluster or a sight in the cluster?
//
//note: if charge is within cluster, most likeyl on which cluster, where most likely to jump off cluster to which site?, dwell time on cluster as single site
//
//
//prob to hop off from a given site in cluster n
//hopoff1 = ProbOnSite1 * Dwell1/sum(dwell) * rate from 1to neighbor n + probonsite1*dwell/sum(dwell)*rate from site1toneigh m
//then normalize hopoff[i]/sum(hopoff)
//
//escape time off cluster
//	1. prob hop off site i time =1/(rate to neigh n + rate to neigh m)
//	escape time = site1hop off time * hop off site 1+...
//
//Tprob_off cluster or escape cluster = tescapsesite1*probhopoffsite1+...
//
//write sample simluation, see test case in matlab
//look at time for maximum efficiency time.h, finde standard deviation and overlap run for multiple interations and find average function time
//
//site pval on and off cluster see matlab code
//
//resolution = arbitray (20)
//
//simluation set up, return time(see matlab), what site it jumped to
//
//
//write test function in library
//
//See matlab code for details
//look to simpligy the equations
//
//pvals internal,
//
//Must id sites, neighbor rates, neighbor ids
//
//also passed visitation freqs and thresh
//
//class data struct sites in cluster (arrays) uses static
// with arrays or linked lists  first elemetn cluster id
//
//time off(tprob off), id which it jumps to
};

int setThresh(int n);
// sets thresh as a static for all functions, not in class cluster


//all other funtions
int potentialCluster(int visitFreq1,int visitFreq2);
	//INPUT: the visitation frequency of site 1 and of site 2
	//OUTPUT: -1 if mal-input or error, 0 if not a cluster, 1 if a cluster

int clusterOrSite(int siteId1, int siteId2);
	//INPUT: the cluster Id of one site, the cluster Id of another site
	//OUTPUT: -1 if mal-input or error, 1 if site to site interaction, 2 if site-cluster, 3 if cluster-cluster

//int neighSiteCluster(site * site1, site * site2, int * neighCluster);
	//INPUT: neighbors ids of sites 1 and 2, an array to write all the neighbors to
	//OUTPUT: -1 if mal-input or error, 1 if successful


// test function and simulation fuction

//int testCluster(void);
//runs all test cases of every function, prints out correct output and what the function prints like such f: (correct) real
//any errors will be reported in the error stream with the funciton call


#endif
