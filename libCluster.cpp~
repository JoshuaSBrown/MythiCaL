#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <memory>

#include "cluster.h"
using namespace std;
/*
#ifdef _E_
#define Err 1
#else
#define Err 0
#endif

static int clusterIdCounter = 0;
static int thresh = 0;

int setThresh(int n){
    if(n>0){
        thresh=n;
        if(Err==1) cout<<"Thresh hold is: "<<thresh<<endl;
        return 1;
    }else{
        if(Err==1) cerr<<"ERROR in setThresh"<<endl;
        return 0;
    }
}

int getThresh(){
	return thresh;
}
/*
cluster::cluster(int            siteId1  //simpler constructer and use an add function
 		 vector<int>    neighId1
		 vector<double> neighRates1
		 int            sizen1
 		 int            visitFreq1
		 int            siteId2
		 vector<int>    neighId2
 		 vector<double> neighRates2
		 int            sizen2
		 int            visitFreq2){
		 
cluster::cluster(){
	clusterId=clusterIdCounter;
	clusterIdCounter++;
	//site  * tmp = new site(siteId1,neighRates1,neighId1,sizen1,visitFreq1);
	//site * tmp2 = new site(siteId2,neighRates2,neighId2,sizen2,visitFreq2);
//	shared_ptr<site> tmp = site(siteId1,vector<double> neighRates,vector<int> neighId1,visitFreq1);
	//create cluster then add
//	sitesInCluster.push_back(tmp);
//	sitesInCluster.push_back(tmp2);	
}

cluster::~cluster(){
}

int cluster::addSite(shared_ptr<site> newSite){ //call potential cluster in add
	//potentialCluster(
	sitesInCluster.push_back(newSite);
	return 1;
}

int cluster::getClusterId(){
	return clusterId;
}

/*
site ** cluster::getSitesInCluster(){
	return siteInCluster;
}

*/

/*
int testCluster(){
	printf("potentialCluster(10,-10,10): (-1) %d\n",potentialCluster(10,-10));
	printf("potentialCluster(10,6,8): (0) %d\n",potentialCluster(10,6,8));
	printf("potentailCluster(10,10,5): (1) %d\n",potentialCluster(10,10,5));
	printf("clusterOrSite(-1,-1): (1) %d\n",clusterOrSite(-1,-1));
	printf("clusterOrSite(-10,2): (-1) %d\n",clusterOrSite(-10,2));
	printf("clusterOrSite(1,-1): (2) %d\n",clusterOrSite(1,-1));
	printf("clusterOrSite(1,1): (3) %d\n",clusterOrSite(1,1));
	return 0;
}



int cluster::printClusterInfo(){
	cout<<"Cluster ID: "<<clusterId<<endl;
	cout<<"Visit Frequency to Cluster: "<< visitFreqCluster<<endl;
	cout<<"Sites in Cluster: "<<endl;
	for(int i = 0; i < sitesInCluster.size(); i++){
		cout<<"Site Id: "<<sitesInCluster[i]->siteId<<"Visit Frequency: "<<sitesInCluster[i]->visitFreq<<endl;
	}
	return 1;
}


int potentialCluster(int visitFreq1, int visitFreq2){
	if(visitFreq1<0||visitFreq2<0||thresh<0){
		fprintf(stderr,"ERROR in potentialCluster\n");//find a way to turn off in preprocessor
		return -1;
	}
	if((visitFreq1 > thresh) && (visitFreq2 > thresh))
		return 1;
	return 0;
}

int clusterOrSite(int clusterId1, int clusterId2){
	if(clusterId1<-1||clusterId2<-1)
		return -1;
	if(clusterId1 == -1 && clusterId2 == -1)
		return 1;
	if(clusterId1 == -1 || clusterId2 == -1)
		return 2;
	return 3;
}

/*
int cluster::neighSiteCluster(site * site1, site * site2, int * neighCluster){
	int i=0, j=0;
	if(site1 ==NULL ||site2==NULL||neighCluster==NULL){
		if(Err) cerr<<"ERROR in neighSiteCluster"<<endl;
		return -1;
	}
	while(site1->neighIds[i]!=-1){ //Use -1 as a terminator
		neighCluster[i]=site1->neighIds[i];
		i++;
	}
	while(site2->neighIds[j]!=-1){
		neighCluster[i+j]=site2->neighIds[j];
		j++;
	}
	return 1;
}
*/
/*
double cluster::dwellTime(){
	site * tmp;
	double sum;
	for(int i =0; i < sitesInCluster.size(); i++){
		tmp=sitesInCluster[i];
		for(int j =0; j < tmp->sizenId; j++){
			sum =+ neighRates[j];
		}
	}
	return (1/sum);
}

double cluster::probHop(site * shipping, site * recieving){ //return an array?
	for(int i = 0; i++; i< SitesInCluster.size()){ //check
		site * tmp = SitesIncluser.pop();
		if(tmp == recieving){
			if(Err) cout<<"Target not conencted to cluster"<<endl;
			return -1;
		}
	}
	for(int i = 0; i++;i<shipping->sizenId){ //combine all for loops
		if(recieving->sId==shipping->nId[i]) int count = i;
	}
	int totalRates;
	for(int i =0; i++; i<shipping->sizenId){
		totalRates+=shipping->nRates[i];
	}
	return shipping->nRates[count]/totalRates;
}

double cluster::hopOffCluster(site * target){
	site * tmp;
	int totalRates=0;
	for(int i =0;i++;i<SitesInCluster.size()){ //figure out pops, pushbacks, and .size() thing
		tmp=SitesInCluster.pop(); 
		for(int j =0; j++; j< sizen){
			totalRates+=nRates[i]; //could calculate at instantation?
		}
	}
	//find way to get target to cluster rate
	return targetToCluster/totalRates
}

double cluster::escapeTime(site * jail){
*/

