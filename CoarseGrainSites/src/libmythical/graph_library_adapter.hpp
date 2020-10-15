
#ifndef MYTHICAL_GRAPH_LIBRARY_ADAPTER_HPP
#define MYTHICAL_GRAPH_LIBRARY_ADAPTER_HPP

#include <iostream>
#include <memory>
#include <string>

#include "site_container.hpp"
#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"

namespace mythical {
  
  // Can only take arguments of type container<unique_ptr<Edge>>
  template<typename T> 
  T convertSitesOutgoingRatesToUniqueWeightedEdges(
      Site_Container site_container, 
      int siteId)
  {
    T container;
    Rate_Map rate_map = site_container.getRates();
    for ( auto neigh : rate_map[siteId] ){
      int neigh_id = neigh.first;
      double * rate = neigh.second;
    
      auto edge_ptr = std::unique_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,*rate));

      container.insert(container.begin(),std::move(edge_ptr));
    }
    return container; 
  }
 
  // Can only take arguments of type container<shared_ptr<Edge>>
  template<typename T> 
  T convertSitesOutgoingRatesToSharedWeightedEdges(
      Site_Container site_container, 
      int siteId)
  {
    T container;
    Rate_Map rate_map = site_container.getRates();
    for ( auto neigh : rate_map[siteId] ){
      int neigh_id = neigh.first;
      double * rate = neigh.second;
    
      auto edge_ptr = std::shared_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,*rate));

      container.insert(container.begin(),std::move(edge_ptr));
    }
    return container; 
  }

  // Same as the above method but for a vector of integers
  template<typename T> 
  T convertSitesOutgoingRatesToSharedWeightedEdges(
      Site_Container site_container, 
      std::vector<int> siteIds)
  {
    T container;
    Rate_Map rate_map = site_container.getRates();
    for(auto siteId : siteIds ){
      for ( auto neigh : rate_map[siteId] ){
        int neigh_id = neigh.first;
        double * rate = neigh.second;

        auto edge_ptr = std::shared_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,*rate));

        container.insert(container.begin(),std::move(edge_ptr));
      }
    }
    return container; 
  }

  template<typename T> 
  T convertSitesOutgoingRatesToTimeSharedWeightedEdges(
      Site_Container site_container, 
      std::vector<int> siteIds)
  {
    T container;
    Rate_Map rate_map = site_container.getRates();
    for(auto siteId : siteIds ){
      for ( auto neigh : rate_map[siteId] ){
        int neigh_id = neigh.first;
        double  time = 1.0/(*neigh.second);
        auto edge_ptr = std::shared_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,time));

        container.insert(container.begin(),std::move(edge_ptr));
      }
    }
    return container; 
  }

  std::unordered_map<int,std::shared_ptr<ugly::GraphNode<std::string>>> 
    convertSitesToEmptySharedNodes(std::vector<int> siteIds);

}

#endif // MYTHICAL_GRAPH_LIBRARY_ADAPTER_HPP

