#ifndef MYTHICAL_BASIN_EXPLORER_HPP
#define MYTHICAL_BASIN_EXPLORER_HPP

#include <vector>

#include "site_container.hpp"
#include "cluster_container.hpp"
#include "../../../UGLY/include/ugly/graphvisitor/graphvisitor_largest_known_value.hpp"

namespace mythical {

class BasinExplorer{
  public:
    BasinExplorer() : threshold_(0.95), max_exploration_count_(5) {};
    void setThreshold(double threshold);
    void setMaxExplorationCount(int count);
    std::vector<int> findBasin(Site_Container& sites,Cluster_Container& clusters, int siteId);
  private:
    double threshold_;
    double fastest_rate_;
    double slowest_rate_;
    double current_sites_fastest_rate_; 
    size_t max_exploration_count_;

    bool rateFastEnough_(double rate);
    double getRate_(Site_Container& sites, std::weak_ptr<ugly::Edge> edge, int vertex);
    void updateFastestRate_(double rate);
    void updateSlowestRate_(double rate);

    void addEdges_(
        Site_Container& sites,
        std::vector<std::weak_ptr<ugly::Edge>> edges_weak, 
        int vertex, 
        ugly::GraphVisitorLargestKnownValue & gv_largest_known);
};


}

#endif // MYTHICAL_BASIN_EXPLORER_HPP
