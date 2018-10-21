
#ifndef GRAPH_EDGEDIRECTED_HPP
#define GRAPH_EDGEDIRECTED_HPP

#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <utility>

#include "../../src/libugly/edge/edge.hpp"

namespace ugly {
  // Composed of two integers describing a link
  // between two vertices
  class EdgeDirected : public Edge {
    public:
      EdgeDirected() : Edge() {
        object_type_ = constants::EdgeType::directed;     
      }
      EdgeDirected(int vertex1, int vertex2) : Edge(vertex1,vertex2) {
        object_type_ = constants::EdgeType::directed;     
      }

      static constants::EdgeType getClassType();
    private:
      static const constants::EdgeType class_type_;
  };
}
#endif // GRAPH_EDGEDIRECTED_HPP
