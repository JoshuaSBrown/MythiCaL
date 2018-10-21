
#include <iostream>
#include <vector>
#include <list>
#include <stdexcept>

#include "../../../include/ugly/edge_weighted.hpp"

using namespace std;

namespace ugly {
    
  const constants::EdgeType EdgeWeighted::class_type_ = constants::EdgeType::weighted;

  constants::EdgeType EdgeWeighted::getClassType() { return EdgeWeighted::class_type_; }

  EdgeWeighted& EdgeWeighted::operator=(const EdgeWeighted & edgeweighted){
    vertex1_ = edgeweighted.vertex1_;
    vertex2_ = edgeweighted.vertex2_;
    object_type_ = edgeweighted.object_type_;
    edge_directed_ = edgeweighted.edge_directed_;
    weight_ = edgeweighted.weight_;
    return *this;
  }
}
