
#include <chrono>

#include "kmc_topology_feature.hpp"

using namespace std;

namespace kmccoarsegrain {

  void occupyTopology_(KMC_TopologyFeature* feature) {
    ++(feature->occupied_);
    ++(feature->total_visit_freq_);
  }

  void occupyTopology_(KMC_TopologyFeature* feature,const int&) {
    ++(feature->occupied_);
    ++(feature->total_visit_freq_);
  }

  void vacateTopology_(KMC_TopologyFeature* feature){
    --(feature->occupied_);
  }

  void vacateTopology_(KMC_TopologyFeature* feature,const int&){
    --(feature->occupied_);
  }

  bool isOccupiedTopology_(const KMC_TopologyFeature* feature){
    return feature->occupied_>0;
  }

  bool isOccupiedTopology_(const KMC_TopologyFeature* feature,const int&){
    return feature->occupied_>0;
  }

  void removeWalker_(KMC_TopologyFeature*,const int&){
    return;
  }

  KMC_TopologyFeature::KMC_TopologyFeature(){
    auto seed = chrono::system_clock::now().time_since_epoch().count();
    random_engine_ = mt19937(seed);
    random_distribution_ = uniform_real_distribution<double>(0.0, 1.0);
    occupied_ = 0;
    escape_time_constant_ = 0.0;
    total_visit_freq_ = 0;
  
    occupy_ptr_ = &occupyTopology_;
    occupy_siteId_ptr_ = &occupyTopology_;

    vacate_ptr_ = &vacateTopology_;
    vacate_siteId_ptr_ = &vacateTopology_;

    isOccupied_ptr_ = &isOccupiedTopology_;
    isOccupied_siteId_ptr_ = &isOccupiedTopology_;

    remove_ptr_ = &removeWalker_;
  }

  void KMC_TopologyFeature::setRandomSeed(const unsigned long seed){
    random_engine_ = mt19937(seed);
  }

  double KMC_TopologyFeature::getDwellTime(const int & ) {
    double number = random_distribution_(random_engine_);
    return (-1.0)*log(number) * escape_time_constant_;
  }

}

