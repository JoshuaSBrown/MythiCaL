
#ifndef UGLY_GRAPHVISITOR_HPP
#define UGLY_GRAPHVISITOR_HPP

#include <vector>
#include <list>
#include <set>
#include <map>
#include <memory>
#include <stdexcept>
#include <functional>

#include "../edge/edge.hpp"
#include "../../../include/ugly/constants.hpp"
#include "../../../include/ugly/edge_directed.hpp"
#include "../../../include/ugly/edge_weighted.hpp"
#include "../../../include/ugly/edge_undirected.hpp"
namespace ugly {

  class GraphVisitor{
    public:
      GraphVisitor() : starting_vertex_set_(false) {
        allowed_edge_types_ = {constants::EdgeType::edge };
      }
      void setStartingVertex(int vertex);

      bool exploredVertex(int vertex);

      template<typename T>  
      std::weak_ptr<T> getNextEdge();

      void exploreEdge(std::weak_ptr<Edge> edge);

      void addEdge(std::weak_ptr<Edge> edge);

      virtual void addEdges(std::vector<std::weak_ptr<Edge>> edges);

      int chooseSourceVertex(std::weak_ptr<Edge> edge) const;
      int chooseTerminalVertex(std::weak_ptr<Edge> edge) const;

      int getUnexploredVertex(std::weak_ptr<Edge> edge) const;
      int getExploredVertex(std::weak_ptr<Edge> edge) const;
      /**
       * \brief Determine if an edge can be added
       *
       * If the edge is directed than the first vertex of the edge corresponds
       * to the source and the second to the drain: souce -> drain
       * A directed edge can only be added if a vertex has already been explroed
       * that is in the drain. For an undirected edge at least one of the 
       * vertices must have been explored. 
       **/
      bool edgeCanBeAdded(std::weak_ptr<Edge> edge);

      /**
       * \brief Determine if both vertices in the edge have been explored
       **/
      bool verticesHaveBeenExplored(std::weak_ptr<Edge> edge) const;
      bool allEdgesExplored() const { return edges_to_explore_.size()==0 ;}
    protected:
      // First int is the vertex, the double is the distance
      std::map<int,double> explored_vertices_;
      std::set<std::weak_ptr<Edge>> explored_edges_;
      // If I do not use a reference wrapper here than I will be unable to take
      // advantage of polymorphism as the list will simply use the Edge class
      std::list<std::weak_ptr<Edge>> edges_to_explore_;
      std::list<constants::EdgeType> allowed_edge_types_;
      
      std::map<constants::EdgeType,std::vector<constants::EdgeType>> allowed_conversions_;

      int starting_vertex_;
      bool starting_vertex_set_;
      bool canAddEdge_(std::weak_ptr<Edge> edge) const;
      virtual void addEdge_(std::weak_ptr<Edge> edge);

      virtual void exploreEdge_(std::weak_ptr<Edge> edge);
      
      bool potentialEdgeKnown_(std::weak_ptr<Edge> edge) const;
      
      bool edgeTypeAllowed_(std::weak_ptr<Edge> edge) const;
      

      template<typename T>
      std::vector<std::weak_ptr<T>> getExploredEdges_();

      virtual std::weak_ptr<Edge> getNextEdge_();
  };

  template<typename T>
  std::weak_ptr<T> GraphVisitor::getNextEdge() {
    auto edge_ptr = getNextEdge_();
    if(auto edge = edge_ptr.lock()){
      if(T::getClassType()==constants::EdgeType::edge){
        return std::static_pointer_cast<T>(edge);
      }
      if(T::getClassType()==edge->getEdgeType()){
        return std::static_pointer_cast<T>(edge);
      }else{
        for( auto type : allowed_conversions_[edge->getEdgeType()] ){
          if(T::getClassType()==type){
            return std::static_pointer_cast<T>(edge);
          }
        }
      }
    }
    throw std::runtime_error("Error cannot retrive edge of the type specified.");
  }
}

#endif // UGLY_GRAPHVISITOR_HPP
