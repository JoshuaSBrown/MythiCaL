
#ifndef UGLY_GRAPHVISITORDEPTHVIRST_HPP
#define UGLY_GRAPHVISITORDEPTHVIRST_HPP

#include "graphvisitor.hpp"

namespace ugly {

  class EdgeWeighted;

  class GraphVisitorDepthFirst : public GraphVisitor {
    public:
      GraphVisitorDepthFirst();
      void addEdges(std::vector<std::weak_ptr<Edge>> edges);
      double getDistanceOfVertex(int vertex);
    private:
      void exploreEdge_(std::weak_ptr<Edge> edge_ptr);
      void addEdge_(std::weak_ptr<Edge> edge_ptr);
      std::weak_ptr<Edge> getNextEdge_();
      std::weak_ptr<Edge> getEdgeShortestDistance_();
      double getDistanceFromStartingVertexToEdge_(std::weak_ptr<Edge> edge);

  };

}

#endif // UGLY_GRAPHVISITORDEPTHFIRST_HPP
