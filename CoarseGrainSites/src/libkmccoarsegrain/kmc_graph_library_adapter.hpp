
#ifndef KMCCOARSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP
#define KMCCOARSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP

#include <iostream>
#include <memory>
#include <string>

#include "kmc_site_container.hpp"
#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"

namespace kmccoarsegrain {
  
  // Can only take arguments of type container<unique_ptr<Edge>>
  template<typename T> 
  T convertSitesOutgoingRatesToUniqueWeightedEdges(
//      KMC_Site_Container & site_container,  
      std::unordered_map<int,double> & neigh_and_rates, 
      const int & siteId)
  {
    T container;
    //Rate_Map rate_map = site_container.getRates();
    //for ( auto& neigh : *rate_map[siteId] ){
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
  T convertSitesOutgoingRatesToSharedWeightedEdges(
      //KMC_Site_Container & site_container, 
      std::unordered_map<int,double> & neigh_and_rates, 
      const int & siteId)
  {
    T container; 
    std::cout << "Getting rates " << std::endl;
//    Rate_Map rate_map = site_container.getRates();
    std::cout << "got rate for site " << siteId << std::endl;
 //   for ( std::pair<int,double>  neigh : *rate_map[siteId] ){
    for ( auto& neigh : neigh_and_rates ){
      std::cout << "site id " << siteId << " neigh id " << neigh.first << std::endl;
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
      //KMC_Site_Container & site_container, 
      KMC_Dynamic_Topology & topology,
      std::vector<int> siteIds)
  {
    T container;
    std::cout << "2 Getting rates " << std::endl;
    //Rate_Map rate_map = site_container.getRates();
    std::cout << "got rate " << std::endl;
    for(auto siteId : siteIds ){
     // for ( auto & neigh : *rate_map[siteId] ){
      for ( auto& neigh : topology.getRates(siteId)){
      std::cout << "site id " << siteId << " neigh id " << neigh.first << std::endl;
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
      //KMC_Site_Container & site_container, 
      KMC_Dynamic_Topology & topology,
      std::vector<int> siteIds)
  {
    T container;
    //Rate_Map rate_map = site_container.getRates();
    for(auto siteId : siteIds ){
     // for ( auto & neigh : *rate_map[siteId] ){
    for ( auto& neigh : topology.getRates(siteId) ){
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

