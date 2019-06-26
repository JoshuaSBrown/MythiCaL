
#include <cassert>
#include <iostream>

#include "../../libkmccoarsegrain/kmc_cluster_container.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

  cout << "Testing: Constructor" << endl;
  {
    KMC_Cluster_Container cluster_container;
  }

  cout << "Testing: addKMC_Cluster" << endl;
  {
    KMC_Cluster_Container cluster_container;
    
    KMC_Cluster cluster;
    cluster.setId(1);

    cluster_container.addKMC_Cluster(cluster);

    bool throw_error = false;
    try {
      cluster_container.addKMC_Cluster(cluster);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }
/*
  cout << "Testing: addKMC_Clusters" << endl;
  {
    KMC_Cluster cluster;
    cluster.setId(1);

    vector<KMC_Cluster> clusters;
    clusters.push_back(cluster);

    KMC_Cluster_Container cluster_container;
    cluster_container.addKMC_Clusters(clusters);

    // Attempt to add the same cluster twice
    clusters.push_back(cluster);
    KMC_Cluster_Container cluster_container2;
    
    bool throw_error = false;
    try{
      cluster_container2.addKMC_Clusters(clusters);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }*/

  cout << "Testing: getKMC_Cluster" << endl;
  {
    KMC_Cluster cluster;
    cluster.setId(1);

    KMC_Cluster_Container cluster_container;

    cluster_container.addKMC_Cluster(cluster);

    auto cluster2 = cluster_container.getKMC_Cluster(1);
    assert(cluster2.getId()==1);

    // Try to grab a cluster that is not stored in the container
    bool throw_error = false;
    try { 
      cluster_container.getKMC_Cluster(0);
    }catch(...) {
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: size" << endl;
  {
    KMC_Cluster cluster;
    KMC_Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    KMC_Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addKMC_Cluster(cluster);
    cluster_container.addKMC_Cluster(cluster2);
    assert(cluster_container.size()==2);
  }

  cout << "Testing: exist" <<endl;
  {
    KMC_Cluster cluster;
    KMC_Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    KMC_Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addKMC_Cluster(cluster);
    cluster_container.addKMC_Cluster(cluster2);
    assert(cluster_container.size()==2);

    assert(cluster_container.exist(1));
    assert(cluster_container.exist(2));
    assert(cluster_container.exist(0)==false);
  }

  cout << "Testing: isOccupied" << endl;
  {
    KMC_Cluster cluster;
    KMC_Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    KMC_Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addKMC_Cluster(cluster);
    cluster_container.addKMC_Cluster(cluster2);
    assert(cluster_container.isOccupied(1)==false); 
    assert(cluster_container.isOccupied(2)==false); 
  }

  cout << "Testing: occupy" << endl;
  {
    KMC_Cluster cluster;
    KMC_Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    KMC_Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addKMC_Cluster(cluster);
    cluster_container.addKMC_Cluster(cluster2);
    assert(cluster_container.isOccupied(1)==false); 
    cluster_container.occupy(1); 
    assert(cluster_container.isOccupied(1));
  }

  cout << "Testing: vacate" << endl;
  {
    KMC_Cluster cluster;
    KMC_Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    KMC_Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addKMC_Cluster(cluster);
    cluster_container.addKMC_Cluster(cluster2);
    assert(cluster_container.isOccupied(1)==false); 
    cluster_container.occupy(1); 
    assert(cluster_container.isOccupied(1));
    cluster_container.vacate(1);
    assert(cluster_container.isOccupied(1)==false); 
  }
  return 0;
}
