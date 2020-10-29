#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <list>
#include <vector>
#include <set>

#include "../../libmythical/graph_library_adapter.hpp"

using namespace std;
using namespace ugly;
using namespace mythical;

TEST_CASE("Testing: graph library adapter","[unit]"){

  cout << "Testing: convertSitesOutgoingRatesToUniqueWeightedEdges" << endl;
  {
    Site site1;
    Site site2;
    Site site3;

    site1.setId(1);
    site2.setId(2);
    site3.setId(3);

    double rate1 = 1.0;
    double rate2 = 2.0;
    double rate3 = 2.0;
    double rate4 = 1.0;

    unordered_map<int, double> neigh_rates_site1;
    neigh_rates_site1[2] = rate1;
    site1.setRatesToNeighbors(neigh_rates_site1);

    unordered_map<int, double> neigh_rates_site2;
    neigh_rates_site2[1] = rate2;
    neigh_rates_site2[3] = rate3;
    site2.setRatesToNeighbors(neigh_rates_site2);
     
    unordered_map<int, double> neigh_rates_site3;
    neigh_rates_site3[2] = rate4;
    site3.setRatesToNeighbors(neigh_rates_site3);

    Site_Container site_container;
    site_container.addSite(site1);
    site_container.addSite(site2);
    site_container.addSite(site3);
   
    cout << "Test with: list" << endl;
    { 
      auto edges = convertSitesOutgoingRatesToUniqueWeightedEdges<list<unique_ptr<Edge>>>(site_container,2);

      assert(edges.size()==2);

      bool found_edge2_1 = false; 
      bool found_edge2_3 = false; 
      for(auto& edge : edges ){
        if(edge->getVertex1()==2){
          if(edge->getVertex2()==1) found_edge2_1 = true;
          if(edge->getVertex2()==3) found_edge2_3 = true;
        }
      }
      assert(found_edge2_1);
      assert(found_edge2_3);
    }

    cout << "Test with: vector" << endl;
    { 
      auto edges = convertSitesOutgoingRatesToUniqueWeightedEdges<vector<unique_ptr<Edge>>>(site_container,2);

      assert(edges.size()==2);

      bool found_edge2_1 = false; 
      bool found_edge2_3 = false; 
      for(auto& edge : edges ){
        if(edge->getVertex1()==2){
          if(edge->getVertex2()==1) found_edge2_1 = true;
          if(edge->getVertex2()==3) found_edge2_3 = true;
        }
      }
      assert(found_edge2_1);
      assert(found_edge2_3);
    }

    cout << "Test with: set" << endl;
    { 
      auto edges = convertSitesOutgoingRatesToUniqueWeightedEdges<set<unique_ptr<Edge>>>(site_container,2);

      assert(edges.size()==2);

      bool found_edge2_1 = false; 
      bool found_edge2_3 = false; 
      for(auto& edge : edges ){
        if(edge->getVertex1()==2){
          if(edge->getVertex2()==1) found_edge2_1 = true;
          if(edge->getVertex2()==3) found_edge2_3 = true;
        }
      }
      assert(found_edge2_1);
      assert(found_edge2_3);
    }
  }
}
