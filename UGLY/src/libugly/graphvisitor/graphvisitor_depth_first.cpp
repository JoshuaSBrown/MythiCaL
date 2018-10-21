
#include <memory>
#include "graphvisitor_depth_first.hpp"
#include "../../../include/ugly/edge_weighted.hpp"

using namespace std;

namespace ugly {

  GraphVisitorDepthFirst::GraphVisitorDepthFirst(){
    allowed_edge_types_.clear();
    allowed_edge_types_.push_back(constants::EdgeType::weighted);
    allowed_edge_types_.push_back(constants::EdgeType::directed_weighted);
    allowed_conversions_[constants::EdgeType::weighted].push_back(constants::EdgeType::directed_weighted);
    allowed_conversions_[constants::EdgeType::directed_weighted].push_back(constants::EdgeType::weighted);
  }

  void GraphVisitorDepthFirst::addEdges(vector<weak_ptr<Edge>> edges){
    for( auto edge_ptr : edges){
      addEdge(edge_ptr);
    }
  }

  double GraphVisitorDepthFirst::getDistanceOfVertex(int vertex){
    if(explored_vertices_.count(vertex)){
      return explored_vertices_[vertex];
    }
    throw invalid_argument("Cannot get distance to vertex as the vertex has not"
        " yet been explored.");
  }

  /****************************************************************************
   * Private Internal Methods
   ****************************************************************************/

  void GraphVisitorDepthFirst::addEdge_(weak_ptr<Edge> ){
    return;
  }

  void GraphVisitorDepthFirst::exploreEdge_(weak_ptr<Edge> edge_ptr){
    auto distance = getDistanceFromStartingVertexToEdge_(edge_ptr);
    if(!verticesHaveBeenExplored(edge_ptr)){
      explored_vertices_[getUnexploredVertex(edge_ptr)]=distance;
    }
  }

  weak_ptr<Edge> GraphVisitorDepthFirst::getEdgeShortestDistance_(){
    double distance;
    bool edge_uninitialized = true;
    auto edge_shortest_distance_away = edges_to_explore_.begin();

    for(auto edge_it=edges_to_explore_.begin();
        edge_it!=edges_to_explore_.end();
        ++edge_it ){

        auto edge_distance = getDistanceFromStartingVertexToEdge_(*edge_it);
        if(edge_uninitialized){
          distance = edge_distance;
          edge_shortest_distance_away = edge_it;
          edge_uninitialized=false;
        }else if(edge_distance<distance){
          distance = edge_distance;
          edge_shortest_distance_away = edge_it;
        }
    }
    if(edge_uninitialized){
      throw runtime_error("Cannot grab edge shortest distance away because "
          "there are no more edges to be explored");
    }
   
    return *edge_shortest_distance_away;
  }

  weak_ptr<Edge> GraphVisitorDepthFirst::getNextEdge_(){
    return GraphVisitorDepthFirst::getEdgeShortestDistance_();
  }

  double 
    GraphVisitorDepthFirst::getDistanceFromStartingVertexToEdge_(
        weak_ptr<Edge> edge_ptr){

      double edge_distance = 0.0;
      if(verticesHaveBeenExplored(edge_ptr)){
        if(auto edge = edge_ptr.lock()){
          auto edge_distance1 = explored_vertices_[edge->getVertex1()];
          auto edge_distance2 = explored_vertices_[edge->getVertex2()];
          if(edge_distance1<edge_distance2){
            edge_distance +=edge_distance1;
          }else{
            edge_distance +=edge_distance2;
          }
        }else{
          throw runtime_error("Cannot lock edge to get vertices of the edge");
        }
      }else{
        auto exploredVertex = getExploredVertex(edge_ptr);
        edge_distance += explored_vertices_[exploredVertex];
        if(!edgeTypeAllowed_(edge_ptr)){
          throw invalid_argument("Edge type is not allowed in the depth first graph"
              " visitor.");
        }
      }

      if( auto edge = edge_ptr.lock()){
        auto edge_weighted_ptr = static_pointer_cast<EdgeWeighted>(edge); 
        edge_distance += edge_weighted_ptr->getWeight();
      }else{
        throw runtime_error("Cannot lock edge to get weight of the edge");
      }
    return edge_distance;
  }
}

