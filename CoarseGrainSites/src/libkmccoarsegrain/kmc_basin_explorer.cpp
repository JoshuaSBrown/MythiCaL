
#include <unordered_set>

#include "kmc_basin_explorer.hpp"
#include "kmc_graph_library_adapter.hpp"
#include "../../../UGLY/include/ugly/graphvisitor/graphvisitor_largest_known_value.hpp"
#include "../../include/kmccoarsegrain/kmc_coarsegrainsystem.hpp"
using namespace ugly;
using namespace std;

namespace kmccoarsegrain {

  typedef unordered_set<shared_ptr<Edge>> shared_edge_set;
  typedef vector<weak_ptr<Edge>> weak_edge_vec;

/*  vector<int> BasinExplorer::findBasin(
      KMC_Site_Container& ,
      KMC_Cluster_Container& ,
      KMC_CoarseGrainSystem & sys,
      int siteId){*/
  vector<int> BasinExplorer::findBasin(
      KMC_Dynamic_Topology & topology, 
      int siteId){
   
   cout << "1 siteid " << siteId << endl; 
    auto edges_store = 
      convertASitesOutgoingRatesToSharedWeightedEdges<shared_edge_set>( topology.getSiteRates(siteId), siteId);

    weak_edge_vec edges_weak(edges_store.begin(),edges_store.end());

   cout << "2 " << endl; 
    GraphVisitorLargestKnownValue gv_largest_known;
    gv_largest_known.setStartingVertex(siteId);

    // Create the site if it does not exist
    if(topology.siteExist(siteId)==false){
      //sys.topology_features_func_[siteId].feature(&sys,siteId);
      topology.features[siteId].feature(topology,siteId);
    }
    //fastest_rate_ = sys.sites_->getFastestRateOffSite(siteId);
    fastest_rate_ = topology.getFastestRateOffSite(siteId);
   cout << "3 " << endl; 

    //if(sys.sites_->partOfCluster(siteId)){
    if(topology.partOfCluster(siteId)){
//      auto & cluster = sys.clusters_->getKMC_Cluster(sys.sites_->getClusterIdOfSite(siteId));
      int clusterId = topology.getClusterIdOfSite(siteId);
      slowest_rate_ = topology.getFastestRateOffCluster(clusterId);
//      slowest_rate_ = cluster.getFastestRateOffCluster();
    }else{
      slowest_rate_ = fastest_rate_;
    }
   cout << "4 " << endl; 
    current_sites_fastest_rate_ = fastest_rate_;

//    addEdges_(*sys.sites_,edges_weak,siteId,gv_largest_known);
    addEdges_(topology,edges_weak,siteId,gv_largest_known);

   cout << "5 " << endl; 
    int exploration_count = 1;

    while(gv_largest_known.allEdgesExplored()==false){
      weak_ptr<Edge> next_edge = gv_largest_known.getNextEdge<Edge>();

   cout << "6 " << endl; 
      auto next_vertex = gv_largest_known.chooseTerminalVertex(next_edge);
      //if(sys.sites_->exist(next_vertex)==false){
      if(topology.siteExist(next_vertex)==false){
//        sys.topology_features_func_[next_vertex].feature(&sys,next_vertex);
        topology.features[next_vertex].feature(topology,next_vertex);
        
      }
   cout << "7 " << endl; 
      gv_largest_known.exploreEdge(next_edge);
      auto edges_tmp = convertASitesOutgoingRatesToSharedWeightedEdges<shared_edge_set>( topology.getSiteRates(next_vertex), next_vertex);
      //auto edges_tmp = convertSitesOutgoingRatesToSharedWeightedEdges<shared_edge_set>( *sys.sites_, next_vertex);
      
      weak_edge_vec edges_weak_tmp(edges_tmp.begin(),edges_tmp.end());
      edges_store.insert(edges_tmp.begin(),edges_tmp.end());

   cout << "8 " << endl; 
      //addEdges_(*sys.sites_,edges_weak_tmp,next_vertex,gv_largest_known);
      addEdges_(topology,edges_weak_tmp,next_vertex,gv_largest_known);
 
      ++exploration_count; 
      if(gv_largest_known.countExploredVertices()>max_exploration_count_){
        vector<int> empty_vec;
        return empty_vec;
      }
   cout << "9 " << endl; 
    }
   cout << "10 " << endl; 
    return gv_largest_known.getExploredVertices();

  }

  void BasinExplorer::addEdges_(

      KMC_Dynamic_Topology& topology,
      weak_edge_vec edges_weak,
      int vertex,
      GraphVisitorLargestKnownValue & gv_largest_known){

    for(auto edge : edges_weak){
      if(gv_largest_known.edgeCanBeAdded(edge)){
        double rate = getRate_(topology,edge, vertex);
        updateFastestRate_(rate);

       // current_sites_fastest_rate_ = sites.getFastestRateOffSite(vertex);
        current_sites_fastest_rate_ = topology.getFastestRateOffSite(vertex);
        // The problem with this is it is updating as it accessing nodes.

        if(rateFastEnough_(rate)){
          
          gv_largest_known.addEdge(edge);
          // Update the slowest rate
          updateSlowestRate_(rate);
        }
      }
    }

  }

  void BasinExplorer::setThreshold(double threshold){
    threshold_ = threshold;
  }

  void BasinExplorer::setMaxExplorationCount(int count){
    max_exploration_count_ = count;
  }

  void BasinExplorer::updateFastestRate_(double rate){
    if(rate>fastest_rate_){
      fastest_rate_ = rate;
    }
  }

  void BasinExplorer::updateSlowestRate_(double rate){
    if(rate<slowest_rate_){
      slowest_rate_ = rate;
    }
  }

  bool BasinExplorer::rateFastEnough_(double rate ){
    double coef1 = rate/fastest_rate_;
    double coef2 = (rate-current_sites_fastest_rate_)/current_sites_fastest_rate_;
    double ratio = coef1*coef2;
    if(ratio>threshold_ || rate >= slowest_rate_*.9999) return true;
    return false;
  }

  double BasinExplorer::getRate_(KMC_Dynamic_Topology & topology,std::weak_ptr<Edge> edge, int vertex){
    double rate;
    if(auto ed = edge.lock()){
      //rate = sites.getRateToNeighborOfSite(vertex,ed->getOtherVertex(vertex));
      rate = topology.getRate(vertex,ed->getOtherVertex(vertex));
    }else{
      throw runtime_error("Edge is no longer accesible in memory.");
    }
    return rate;
  }

}

