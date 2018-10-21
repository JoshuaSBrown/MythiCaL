#ifndef UGLY_EDGE_HPP
#define UGLY_EDGE_HPP

#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <utility>
#include <functional>

#include "../../../include/ugly/constants.hpp"

namespace ugly {
  // Composed of two integers describing a link
  // between two vertices
  class Edge{
    protected:
      bool edge_directed_;
      constants::EdgeType object_type_;
      int vertex1_;
      int vertex2_;
    public:
      Edge() :
        edge_directed_(true),
        vertex1_(constants::unassigned),
        vertex2_(constants::unassigned){
          object_type_ = constants::EdgeType::edge;
        }

      explicit Edge(bool directed) : 
        edge_directed_(directed),
        vertex1_(constants::unassigned),
        vertex2_(constants::unassigned){
          object_type_ = constants::EdgeType::edge;
        }
  
      Edge(int vertex1, int vertex2) :
        edge_directed_(true), 
        vertex1_(vertex1), 
        vertex2_(vertex2){
          object_type_ = constants::EdgeType::edge;
        }


      Edge(const Edge &edge){
        edge_directed_ = edge.edge_directed_, 
        object_type_ = edge.object_type_;
        vertex1_ = edge.vertex1_;
        vertex2_ = edge.vertex2_;
      }

      virtual ~Edge() {};

      bool directional() const { return edge_directed_; }

      int getOtherVertex(int vertex) const;
      int getVertex1(void) const;
      int getVertex2(void) const;

      int getMax(void) const;
      int getMin(void) const;

      bool containsVertex(int vertex) const;

      Edge& operator=(const Edge &edge);

      bool operator==(const Edge& edge) const;
      bool operator!=(const Edge& edge) const;

      // edge 1 is less than edge 2 if the smallest vertex in edge 1 (vert1)
      // is smaller than the vert1 in edge2. If they are the same than edge1
      // is less than edge2 if the second vertex (vert2) is smaller in edge1
      // than in edge2
      bool operator<(const Edge& edge) const;
      bool operator>(const Edge& edge) const;
      bool operator<=(const Edge& edge) const;
      bool operator>=(const Edge& edge) const;
      friend std::ostream& operator<<(std::ostream& os, const Edge edge); 

      constants::EdgeType getEdgeType() const { return object_type_; }
      static constants::EdgeType getClassType();
    private:
      static const constants::EdgeType class_type_;
  };

}

namespace std{
  template<>
    class hash<ugly::Edge>{
      public:
        size_t operator()(const ugly::Edge &edge) const
        {
          return hash<int>()(edge.getVertex1())^hash<int>()(edge.getVertex2());
        }
    };
}

#endif // UGLY_EDGE_HPP
