#ifndef UGLY_EDGEDIRECTEDWEIGHTED_HPP
#define UGLY_EDGEDIRECTEDWEIGHTED_HPP

#include <iostream>
#include <vector>
#include <list>
#include <cassert>
#include <utility>

#include "../../src/libugly/edge/edge.hpp"

namespace ugly {
// Composed of two integers describing a link
// between two vertices
  class EdgeDirectedWeighted : public Edge {
    private:
      static const constants::EdgeType class_type_;
    protected:
      double weight_;
    public:
      EdgeDirectedWeighted() : weight_(1.0) {
          object_type_ = constants::EdgeType::directed_weighted;
          edge_directed_ = true;
        }

      EdgeDirectedWeighted(int vertex1, int vertex2): 
        Edge(vertex1,vertex2) {
          object_type_ = constants::EdgeType::directed_weighted;
          weight_ = 1.0;
        }

      EdgeDirectedWeighted(int vertex1, int vertex2,double weight) :  
        Edge(vertex1,vertex2){
          object_type_ = constants::EdgeType::directed_weighted;
          weight_ = weight;
        }

      EdgeDirectedWeighted(const EdgeDirectedWeighted &edgedirectedweighted) :
        Edge() {
        vertex1_ = edgedirectedweighted.vertex1_; 
        vertex2_ = edgedirectedweighted.vertex2_;
        object_type_ = edgedirectedweighted.object_type_; 
        weight_ = edgedirectedweighted.weight_ ;
        edge_directed_ = true;
      }

      EdgeDirectedWeighted& operator=(const EdgeDirectedWeighted &EdgeDirectedWeighted);

      void setWeight(double weight){ weight_ = weight; }
      double getWeight() const { return weight_; }

      static constants::EdgeType getClassType();

  };

}

#endif // UGLY_EDGEDIRECTEDWEIGHTED_HPP
