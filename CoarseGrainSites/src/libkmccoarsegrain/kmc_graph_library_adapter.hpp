
#ifndef KMCCOARSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP
#define KMCCOARSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP

#include <iostream>
#include <memory>
#include <string>

#include "kmc_dynamic_topology.hpp"
#include "kmc_site_container.hpp"
#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"

namespace kmccoarsegrain {
  
  // Can only take arguments of type container<unique_ptr<Edge>>
  template<typename T> 
  T convertASitesOutgoingRatesToUniqueWeightedEdges(
      std::unordered_map<int,double> & neigh_and_rates, 
      const int & siteId)
  {
    T container;
    for ( auto& neigh : neigh_and_rates ){
      int neigh_id = neigh.first;
      const double rate = neigh.second;
    
      auto edge_ptr = std::unique_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,rate));

      container.insert(container.begin(),std::move(edge_ptr));
    }
    return container; 
  }
 
  // Can only take arguments of type container<shared_ptr<Edge>>
  template<typename T> 
  T convertASitesOutgoingRatesToSharedWeightedEdges(
      std::unordered_map<int,double> & neigh_and_rates, 
      const int & siteId)
  {
    T container; 
    for ( auto& neigh : neigh_and_rates ){
      int neigh_id = neigh.first;
      const double rate = neigh.second;
    
      auto edge_ptr = std::shared_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,rate));

      container.insert(container.begin(),std::move(edge_ptr));
    }
    return container; 
  }

  // Same as the above method but for a vector of integers
  template<typename T> 
  T convertSitesOutgoingRatesToSharedWeightedEdges(
      KMC_Dynamic_Topology & topology,
      std::vector<int> siteIds)
  {
    T container;
    for(auto siteId : siteIds ){
      for ( auto& neigh : topology.getSiteRates(siteId)){
        int neigh_id = neigh.first;
        const double  rate = neigh.second;

        auto edge_ptr = std::shared_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,rate));

        container.insert(container.begin(),std::move(edge_ptr));
      }
    }
    return container; 
  }

  template<typename T> 
  T convertSitesOutgoingRatesToTimeSharedWeightedEdges(
      KMC_Dynamic_Topology & topology,
      std::vector<int> siteIds)
  {
    T container;
    for(auto siteId : siteIds ){
      for ( auto& neigh : topology.getSiteRates(siteId) ){
        int neigh_id = neigh.first;
        double  time = 1.0/(neigh.second);
        auto edge_ptr = std::shared_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,time));

        container.insert(container.begin(),std::move(edge_ptr));
      }
    }
    return container; 
  }

  std::unordered_map<int,std::shared_ptr<ugly::GraphNode<std::string>>> 
    convertSitesToEmptySharedNodes(std::vector<int> siteIds);


}

#endif // KMCCOARSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP

