
#ifndef KMCCOURSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP
#define KMCCOURSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP

#include <iostream>
#include <memory>

#include "kmc_site_container.hpp"
#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"

namespace kmccoursegrain {
  
/*
  sitesToEdges(KMC_Site_Container site_container); 
  
  sitesToEdges(KMC_Site_Container site_container, vector<int> siteIds); 
*/

  // Can only take arguments of type container<unique_ptr<Edge>>
  template<typename T> 
  T convertSitesOutgoingRatesToWeightedEdges(KMC_Site_Container site_container, int siteId){
    T container;
    Rate_Map rate_map = site_container.getRates();
    for ( auto neigh : rate_map[siteId] ){
      int neigh_id = neigh.first;
      double * rate = neigh.second;
    
      ugly::EdgeDirectedWeighted edgew(siteId,neigh_id,*rate);
      auto edge_ptr = std::unique_ptr<ugly::EdgeDirectedWeighted>(new ugly::EdgeDirectedWeighted(siteId,neigh_id,*rate));

      container.insert(container.begin(),std::move(edge_ptr));
    }
    return container; 
  }
  
}

#endif // KMCCOURSEGRAIN_KMC_GRAPH_LIBRARY_ADAPTER_HPP

