
#include <iostream>
#include <vector>
#include <list>
#include <stdexcept>

#include "../../../include/ugly/constants.hpp"
#include "edge.hpp"

using namespace std;

namespace ugly {

  const constants::EdgeType Edge::class_type_ = constants::EdgeType::edge;

  constants::EdgeType Edge::getClassType() { return Edge::class_type_; }

  int Edge::getOtherVertex(int vertex) const{
    if(vertex==vertex1_) {
      if(vertex2_ == constants::unassigned){
        throw runtime_error("vertex2_ has not been assigned a value"); 
      }
      return vertex2_;
    }else if(vertex==vertex2_){
      if(vertex1_ == constants::unassigned){
        throw runtime_error("vertex1_ has not been assigned a value"); 
      }
      return vertex1_;
    }
    throw runtime_error("Unable to determine the other vertex as the provided "
        "vertex is not stored in the edge.");
  }

  int Edge::getVertex1() const {
    if(vertex1_==constants::unassigned){
      throw runtime_error("Cannot get vertex 1 as it has not been assigned.");
    }
    return vertex1_;
  }

  int Edge::getVertex2() const {
    if(vertex2_==constants::unassigned){
      throw runtime_error("Cannot get vertex 2 as it has not been assigned.");
    }
    return vertex2_;
  }

  int Edge::getMin() const {
    if(vertex1_==constants::unassigned || vertex2_==constants::unassigned){
      throw runtime_error("Cannot get min vertex as at least one vertex in the "
         "edge has not been assigned.");
    }
    return min(vertex1_,vertex2_);
  }

  int Edge::getMax() const {
    if(vertex1_==constants::unassigned || vertex2_==constants::unassigned){
      throw runtime_error("Cannot get max vertex as at least one vertex in the "
         "edge has not been assigned.");
    }
    return max(vertex1_,vertex2_);
  }

  bool Edge::containsVertex(int vertex) const {
    return (vertex1_==vertex || vertex2_==vertex);
  }

  Edge& Edge::operator=(const Edge &edge){
    vertex1_ = edge.vertex1_;
    vertex2_ = edge.vertex2_;
    object_type_ = edge.object_type_;
    edge_directed_ = edge.edge_directed_;
    return *this;
  }

  bool Edge::operator==(const Edge & edge) const{

    if(vertex1_==edge.vertex1_ && vertex2_==edge.vertex2_) return true; 
    if(vertex2_==edge.vertex1_ && vertex2_==edge.vertex2_) return true; 
    return false;
  }

  bool Edge::operator!=(const Edge& edge) const{
    return !(*this==edge);
  }

  bool Edge::operator<(const Edge& edge) const{
    if(vertex1_<edge.vertex1_) return true;
    if(vertex1_>edge.vertex1_) return false;
    if(vertex2_<edge.vertex2_) return true;
    return false;
  }

  bool Edge::operator<=(const Edge& edge) const{
    return (*this<edge || *this==edge);
  }

  bool Edge::operator>(const Edge& edge) const{
    return !(*this<=edge);
  }


  bool Edge::operator>=(const Edge& edge) const{
    return !(*this<edge);   
  }

  ostream& operator<<(ostream& os, const Edge edge){
    os << "Vertices" << endl;
    os << edge.vertex1_ << " " << edge.vertex2_ << endl; 
    return os;
  }

}
