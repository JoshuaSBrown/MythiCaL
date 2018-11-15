
#include "kmc_graph_library_adapter.hpp"

using namespace std;
using namespace ugly;

namespace kmccoursegrain {

  unordered_map<int,shared_ptr<GraphNode<string>>> 
  convertSitesToEmptySharedNodes(vector<int> siteIds)
  {
    unordered_map<int, shared_ptr<GraphNode<string>>> nodes;
    for (auto siteId : siteIds) {
      if (nodes.count(siteId) == 0) {
        auto graph_node = shared_ptr<GraphNode<string>>(new GraphNode<string>(""));
        nodes[siteId] = graph_node;
      }
    }
    return nodes;
  }

}


