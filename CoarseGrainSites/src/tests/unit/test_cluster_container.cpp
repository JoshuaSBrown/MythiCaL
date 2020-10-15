
#include <cassert>
#include <iostream>

#include "../../libmythical/cluster_container.hpp"

using namespace std;
using namespace mythical;

int main(void){

  cout << "Testing: Constructor" << endl;
  {
    Cluster_Container cluster_container;
  }

  cout << "Testing: addCluster" << endl;
  {
    Cluster_Container cluster_container;
    
    Cluster cluster;
    cluster.setId(1);

    cluster_container.addCluster(cluster);

    bool throw_error = false;
    try {
      cluster_container.addCluster(cluster);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: addClusters" << endl;
  {
    Cluster cluster;
    cluster.setId(1);

    vector<Cluster> clusters;
    clusters.push_back(cluster);

    Cluster_Container cluster_container;
    cluster_container.addClusters(clusters);

    // Attempt to add the same cluster twice
    clusters.push_back(cluster);
    Cluster_Container cluster_container2;
    
    bool throw_error = false;
    try{
      cluster_container2.addClusters(clusters);
    }catch(...){
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: getCluster" << endl;
  {
    Cluster cluster;
    cluster.setId(1);

    Cluster_Container cluster_container;

    cluster_container.addCluster(cluster);

    auto cluster2 = cluster_container.getCluster(1);
    assert(cluster2.getId()==1);

    // Try to grab a cluster that is not stored in the container
    bool throw_error = false;
    try { 
      cluster_container.getCluster(0);
    }catch(...) {
      throw_error = true;
    }
    assert(throw_error);
  }

  cout << "Testing: size" << endl;
  {
    Cluster cluster;
    Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addCluster(cluster);
    cluster_container.addCluster(cluster2);
    assert(cluster_container.size()==2);
  }

  cout << "Testing: exist" <<endl;
  {
    Cluster cluster;
    Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addCluster(cluster);
    cluster_container.addCluster(cluster2);
    assert(cluster_container.size()==2);

    assert(cluster_container.exist(1));
    assert(cluster_container.exist(2));
    assert(cluster_container.exist(0)==false);
  }

  cout << "Testing: isOccupied" << endl;
  {
    Cluster cluster;
    Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addCluster(cluster);
    cluster_container.addCluster(cluster2);
    assert(cluster_container.isOccupied(1)==false); 
    assert(cluster_container.isOccupied(2)==false); 
  }

  cout << "Testing: occupy" << endl;
  {
    Cluster cluster;
    Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addCluster(cluster);
    cluster_container.addCluster(cluster2);
    assert(cluster_container.isOccupied(1)==false); 
    cluster_container.occupy(1); 
    assert(cluster_container.isOccupied(1));
  }

  cout << "Testing: vacate" << endl;
  {
    Cluster cluster;
    Cluster cluster2;

    cluster.setId(1);
    cluster2.setId(2);

    Cluster_Container cluster_container;
    assert(cluster_container.size()==0);

    cluster_container.addCluster(cluster);
    cluster_container.addCluster(cluster2);
    assert(cluster_container.isOccupied(1)==false); 
    cluster_container.occupy(1); 
    assert(cluster_container.isOccupied(1));
    cluster_container.vacate(1);
    assert(cluster_container.isOccupied(1)==false); 
  }
  return 0;
}
