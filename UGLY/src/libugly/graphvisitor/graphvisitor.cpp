#include <algorithm>
#include <stdexcept>
#include <iostream>
#include "../edge/edge.hpp"
#include "graphvisitor.hpp"
//#include "../reference_wrapper_suplement.hpp"
#include "../weak_pointer_supplement.hpp"

using namespace std;

namespace ugly {

  /****************************************************************************
   * External public facing methods
   ****************************************************************************/

  void GraphVisitor::setStartingVertex(int vertex){

    starting_vertex_ = vertex;
    explored_vertices_[vertex] = 0.0;
    starting_vertex_set_ = true;
  }

  bool GraphVisitor::exploredVertex(int vertex){
    return explored_vertices_.count(vertex);
  }

  bool GraphVisitor::edgeCanBeAdded(weak_ptr<Edge> edge_ptr) {
    if(find(edges_to_explore_.begin(),
          edges_to_explore_.end(),
          edge_ptr)!=edges_to_explore_.end()){
      return false;
    }
    if(explored_edges_.count(edge_ptr)){
      return false;
    }
    if(auto edge = edge_ptr.lock()){
      if(explored_vertices_.count(edge->getVertex1())) {
        return true; 
      }
      if(!edge->directional()){
        if(explored_vertices_.count(edge->getVertex2())){ 
          return true; 
        }
      }
    }
    return false;
  }

  void GraphVisitor::exploreEdge(weak_ptr<Edge> edge_ptr){
    if(potentialEdgeKnown_(edge_ptr)==false){
      throw invalid_argument("Cannot explore an edge that has not yet been "
          "added to the visitor. Exploration consists of adding a vertex a "
          "vertex that is being explored if it has not been explored yet, or "
          "simply updating the graphnode. It also consists of removing the edge"
          "from the list of potential eddges to explore.");
    }
    if(explored_edges_.count(edge_ptr)){
      throw invalid_argument("Cannot explore edge that has already been "
          "explored.");
    }
    exploreEdge_(edge_ptr);
    explored_edges_.insert(edge_ptr);
    edges_to_explore_.remove(edge_ptr);
  }

  void GraphVisitor::exploreEdge_(weak_ptr<Edge>){
    throw runtime_error("Cannot call exploreEdge_ from base class must define it"
        " in the derived class.");
  }

  void GraphVisitor::addEdge(weak_ptr<Edge> edge_ptr){

    if(explored_edges_.count(edge_ptr)){
      throw invalid_argument("Cannot add edge it has already been explored");
    }
    if( find(edges_to_explore_.begin(),edges_to_explore_.end(),edge_ptr)!=edges_to_explore_.end()){
      throw invalid_argument("Cannot add the same edge more than twice.");
    }
/*    if(auto edge = edge_ptr.lock()){
      for(auto edge_temp_ptr : edges_to_explore_ ){
        if(auto edge_temp = edge_temp_ptr.lock()){
          if(edge_temp==edge) {
            throw invalid_argument("Cannot add the same edge more than twice.");
          }
        }
      }
    }*/
    if(!starting_vertex_set_){
      throw runtime_error("Cannot add edges before specifying the starting "
          "vertex.");
    }
    if(!canAddEdge_(edge_ptr)){
      throw runtime_error("Cannot add edge as non of the vertices in the edge "
          "have been explored. You can only add edges to the visitor that are "
          "connected to a vertex that has been explored. E.g. Say I set my "
          "starting vertex to 0 and want to add edges 0-1 and 1-2, I cannot do "
          "this because none of the vertices in 1-2 have been explored. If the,"
          " edge is a directed edge the source is always vertex 1 thus in the "
          "above case I could still not add vertex 0-1 as it is directional "
          "from vertex 0 going to vertex 1");
    }
    if(!edgeTypeAllowed_(edge_ptr)){
      throw runtime_error("Cannot add edge of this time to the visitor it is "
          "not allowed.");
    }
    addEdge_(edge_ptr);
    edges_to_explore_.push_back(edge_ptr);
  }

  void GraphVisitor::addEdges(vector<weak_ptr<Edge>>){
    throw runtime_error("Cannot call addEdges from base class must define it in"
        " the derived class.");
  }

  int GraphVisitor::chooseTerminalVertex(weak_ptr<Edge> edge_ptr) const {
    if( auto edge = edge_ptr.lock() ){
      if(explored_vertices_.count(edge->getVertex1())){
        return edge->getVertex2();
      }else if(explored_vertices_.count(edge->getVertex2())){
        return edge->getVertex1();
      }
      return edge->getMax();
    }
    throw invalid_argument("edge is no longer safe to use cannot call lock");
  }

  int GraphVisitor::chooseSourceVertex(weak_ptr<Edge> edge_ptr) const {

    if( auto edge = edge_ptr.lock() ){
      if(explored_vertices_.count(edge->getVertex1())){
        return edge->getVertex1();
      }else if(explored_vertices_.count(edge->getVertex2())){
        return edge->getVertex2();
      }
      return edge->getMin();
    }
    throw invalid_argument("edge is no longer safe to use cannot call lock");
  }

  int GraphVisitor::getUnexploredVertex(weak_ptr<Edge> edge_ptr) const {
    if(verticesHaveBeenExplored(edge_ptr)){
      throw invalid_argument("Cannot find unexplored vertex as both vertices "
        "have been explored.");
    }
    if( auto edge = edge_ptr.lock() ){
      if(explored_vertices_.count(edge->getVertex1())){
        return edge->getVertex2();
      }
      if(explored_vertices_.count(edge->getVertex2())){
        return edge->getVertex1();
      }
    }
    throw invalid_argument("Neither vertices have been explored so you cannot "
        "determine which vertex to return when calling getUnexploredVertex.");
  }

  int GraphVisitor::getExploredVertex(weak_ptr<Edge> edge_ptr) const {
    if(auto edge = edge_ptr.lock()){
      return edge->getOtherVertex(getUnexploredVertex(edge_ptr));
    }
    throw invalid_argument("Cannot convert weak pointer to shared pointer.");
  }

  /****************************************************************************
   * Internal private functions
   ****************************************************************************/

  void GraphVisitor::addEdge_(weak_ptr<Edge>){
    throw runtime_error("Cannot call addEdge from base class must define it in "
        "the derived class.");
  }

  bool GraphVisitor::verticesHaveBeenExplored(weak_ptr<Edge> edge_ptr) const {
    if(auto edge = edge_ptr.lock()){
      if(explored_vertices_.count(edge->getVertex1()) && 
          explored_vertices_.count(edge->getVertex2())){ 
        return true;
      }
    }
    return false;
  }

  bool GraphVisitor::edgeTypeAllowed_(weak_ptr<Edge> edge_ptr) const {
    if( auto edge = edge_ptr.lock()){
      for( auto edge_type : allowed_edge_types_){
        if(edge->getEdgeType()==edge_type){
          return true;
        }
      }
    }
    return false;
  }

  bool GraphVisitor::potentialEdgeKnown_(weak_ptr<Edge> edge_ptr) const {
    auto it = find(edges_to_explore_.begin(),edges_to_explore_.end(),edge_ptr);
    return it!=edges_to_explore_.end();
  }

  bool GraphVisitor::canAddEdge_(weak_ptr<Edge> edge_ptr) const {
    if( auto edge = edge_ptr.lock()){
      if(edge->directional()==true){
        // First int is the source it must be part of the explored vertices
        if(explored_vertices_.count(edge->getVertex1())) return true;
        return false;
      }else{
        if(explored_vertices_.count(edge->getVertex1())) return true;
        if(explored_vertices_.count(edge->getVertex2())) return true;
      }
    }
    return false;
  }

  weak_ptr<Edge> GraphVisitor::getNextEdge_(){
    throw runtime_error("Cannot call getNextEdge from base class must "
        "define in the derived class.");
  }
}
