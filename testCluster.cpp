#include <iostream>

#include "cluster.h"

using namespace std;

#ifdef _T_
#define Var 1
#else
#define Var 0
#endif
/*
int testCluster(cluster * c){
	printf("potentialCluster(10,-10,10): (-1) %d\n",c->potentialCluster(10,-10));
	printf("potentialCluster(10,6,8): (0) %d\n",c->potentialCluster(10,6,8));
	printf("potentailCluster(10,10,5): (1) %d\n",c->potentialCluster(10,10,5));
	printf("clusterOrSite(-1,-1): (1) %d\n",c->clusterOrSite(-1,-1));
	printf("clusterOrSite(-10,2): (-1) %d\n",c->clusterOrSite(-10,2));
	printf("clusterOrSite(1,-1): (2) %d\n",c->clusterOrSite(1,-1));
	printf("clusterOrSite(1,1): (3) %d\n",c->clusterOrSite(1,1));
	return 0;
}
*/
int main(void){
    if(Var==1) cout<<"Yes"<<endl;
    else cout<<"No"<<endl;
    setThresh(10);
    return 0;
}
