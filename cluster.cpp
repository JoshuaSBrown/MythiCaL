#include <stdio.h>
#include <stdlib.h>

#include "cluster.h"

static int clusterIdCounter = 0;
static int thresh = 0;

cluster::cluster(int sitedId1, double * neighId1, int sizenId1, int visitFreq1, int siteId2, double * neighId2, int sizenId2, int visitFreq2){
	clusterId=clusterIdCounter;
	clusterIdCounter++;
	thresh = t;
	if(potentialCluster(visitFreq1,visitFreq2){
			siteInCluster[0]= new site(siteId1,neighId1,sizen1,visitFreq1);
			siteInCluster[1]= new site(sitedId2,neighId2,sizen2,visitFreq2);
			}
}
cluster::~cluster(){
}

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

int cluster::potentialCluster(int visitFreq1, int visitFreq2){
	if(visitFreq1<0||visitFreq2<0||thresh<0){
		fprintf(stderr,"ERROR in potentialCluster\n");//find a way to turn off in preprocessor
		return -1;
	}
	if((visitFreq1 > thresh) && (visitFreq2 > thresh))
		return 1;
	return 0;
}

int cluster::clusterOrSite(int clusterId1, int clusterId2){
	if(clusterId1<-1||clusterId2<-1)
		return -1;
	if(clusterId1 == -1 && clusterId2 == -1)
		return 1;
	if(clusterId1 == -1 || clusterId2 == -1)
		return 2;
	return 3;
}

int cluster::neighSiteCluster(int neighIds1[], int * neighIds2, int * neighCluster){
}
