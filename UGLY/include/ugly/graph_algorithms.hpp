#ifndef UGLY_GRAPHALGORITHMS_HPP
#define UGLY_GRAPHALGORITHMS_HPP

#include <vector>
#include <unordered_map> 
#include <algorithm>
#include <stdexcept>

#include "../../src/libugly/edge/edge.hpp"
#include "../../src/libugly/graphvisitor/graphvisitor_depth_first.hpp"
#include "graph.hpp"
#include "graph_node.hpp"
#include "constants.hpp"

namespace ugly {

  namespace graphalgorithms {

    template<class... Ts> 
      double dijkstraGoingFrom(int start_vertex, int end_vertex, Graph<Ts...>& graph){

        GraphVisitorDepthFirst graphvisitor_depth_first;
        auto edges = graph.getEdgesOriginatingFromVertex(start_vertex);

        graphvisitor_depth_first.setStartingVertex(start_vertex);
        graphvisitor_depth_first.addEdges(edges);

        auto next_edge = graphvisitor_depth_first.getNextEdge<EdgeWeighted>();
        while(graphvisitor_depth_first.allEdgesExplored()==false){

          next_edge = graphvisitor_depth_first.getNextEdge<EdgeWeighted>();
          auto next_vertex = graphvisitor_depth_first.chooseTerminalVertex(next_edge);
          graphvisitor_depth_first.exploreEdge(next_edge);
          edges = graph.getEdgesOriginatingFromVertex(next_vertex);

          for (auto ed : edges ){
            if(graphvisitor_depth_first.edgeCanBeAdded(ed)){
              graphvisitor_depth_first.addEdge(ed);
            }
          }

          if(end_vertex==next_vertex){
            return graphvisitor_depth_first.getDistanceOfVertex(end_vertex);
          }

        }

        throw std::runtime_error("No connection has been found between your two "
            "vertices.");
      }


    template<class... Ts>
      std::map<std::pair<int,int>,double> maxMinimumDistanceBetweenEveryVertex( Graph<Ts...>& graph){

        bool graph_directional = graph.directional();

        std::map<std::pair<int,int>,double> maxMinDistanceOfGraphVertices;
        auto vertices = graph.getVertices();

        for(auto vertices_it1 = vertices.begin();
            vertices_it1!=vertices.end();
            ++vertices_it1){

          int source = constants::unassigned;
          int drain = constants::unassigned;
          double maxdistance = 0.0;

          if(graph_directional){

            for(auto vertices_it2=vertices.begin();
                vertices_it2!=vertices.end();
                ++vertices_it2){

              if(vertices_it2!=vertices_it1){
                auto distance = dijkstraGoingFrom(*vertices_it1,*vertices_it2,graph);
                if(distance>maxdistance) {
                  maxdistance=distance;
                  source = *vertices_it1;
                  drain = *vertices_it2;
                }
              }
            }
          }else{
            auto vertices_it2 = vertices_it1;
            for(++vertices_it2;vertices_it2!=vertices.end();++vertices_it2){

              auto distance = dijkstraGoingFrom(*vertices_it1,*vertices_it2,graph);
              if(distance>maxdistance) {
                maxdistance=distance;
                source = *vertices_it1;
                drain = *vertices_it2;
              }
            }
          }
          if(source!=constants::unassigned){
            maxMinDistanceOfGraphVertices[std::pair<int,int>(source,drain)]=maxdistance;
          }
        }

        return maxMinDistanceOfGraphVertices;  
      }
  }
}
#endif // UGLY_GRAPHALGORITHMS_HPP
