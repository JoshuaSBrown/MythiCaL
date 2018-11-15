
#ifndef KMCCOURSEGRAIN_KMC_BASIN_EXPLORER_HPP
#define KMCCOURSEGRAIN_KMC_BASIN_EXPLORER_HPP

#include <vector>

#include "kmc_site_container.hpp"
#include "../../../UGLY/include/ugly/graphvisitor/graphvisitor_largest_known_value.hpp"

namespace kmccoursegrain {

class BasinExplorer{
  public:
    BasinExplorer() : threshold_(0.95), max_exploration_count_(5) {};
    void setThreshold(double threshold);
    void setMaxExplorationCount(int count);
    std::vector<int> findBasin(KMC_Site_Container& sites,int siteId);
  private:
    double threshold_;
    double fastest_rate_;
    double slowest_rate_;
    double current_sites_fastest_rate_; 
    size_t max_exploration_count_;

    bool rateFastEnough_(double rate);
    double getRate_(KMC_Site_Container& sites, std::weak_ptr<ugly::Edge> edge, int vertex);
    void updateFastestRate_(double rate);
    void updateSlowestRate_(double rate);

    void addEdges_(
        KMC_Site_Container& sites,
        std::vector<std::weak_ptr<ugly::Edge>> edges_weak, 
        int vertex, 
        ugly::GraphVisitorLargestKnownValue & gv_largest_known);
};


}

#endif // KMCCOURSEGRAIN_KMC_BASIN_EXPLORER_HPP
