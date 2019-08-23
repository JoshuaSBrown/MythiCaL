
#include <algorithm>
#include <chrono>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <unordered_set>

#include "../../include/kmccoarsegrain/kmc_coarsegrainsystem.hpp"
#include "../../include/kmccoarsegrain/kmc_constants.hpp"
#include "../../include/kmccoarsegrain/kmc_walker.hpp"

#include "log.hpp"
#include "kmc_basin_explorer.hpp"
#include "kmc_graph_library_adapter.hpp"
#include "kmc_cluster_container.hpp"
#include "kmc_dynamic_topology.hpp"

#include "../../../UGLY/include/ugly/pair_hash.hpp"
#include "../../../UGLY/include/ugly/edge_directed_weighted.hpp"
#include "../../../UGLY/include/ugly/graph.hpp"
#include "../../../UGLY/include/ugly/graph_algorithms.hpp"
#include "../../../UGLY/include/ugly/graph_node.hpp"
#include "../../../UGLY/include/ugly/graphvisitor/graphvisitor_smallest_known_value.hpp"

using namespace std;
using namespace std::chrono;
using namespace ugly;
using namespace ugly::graphalgorithms;

namespace kmccoarsegrain {
  /****************************************************************************
   * Private Internal Function Declarations
   ****************************************************************************/

  bool compare(const pair<int,double> &x, const pair<int,double> &y){
    return x.second>y.second;
  }

  unordered_map<int, shared_ptr<GraphNode<string>>> createNode_(int siteIds);

  size_t countUniqueClusters(const unordered_map<int,int> & sites_and_clusters);
  int getFavoredClusterId(unordered_map<int,int> sites_and_clusters);

  /****************************************************************************
   * Public Facing Functions
   ****************************************************************************/

  KMC_CoarseGrainSystem::KMC_CoarseGrainSystem() :
    performance_ratio_(1.40),
    time_resolution_set_(false),
    minimum_coarse_graining_resolution_(2),
    iteration_(0),
    iteration_threshold_(1000),
    iteration_threshold_min_(1000){
      topology_ = unique_ptr<KMC_Dynamic_Topology>(new KMC_Dynamic_Topology);
    }

  KMC_CoarseGrainSystem::~KMC_CoarseGrainSystem(){
  }

  double KMC_CoarseGrainSystem::getTimeResolution() { 
    if(!time_resolution_set_){
      throw runtime_error("Cannot get the time resolution as it has not yet "
          "been set.");
    }
    return time_resolution_; 
  }

  void KMC_CoarseGrainSystem::setTimeResolution(double time_resolution){
    if(time_resolution<=0.0){
      throw invalid_argument("The time resolution must be a positive value.");
    }
    time_resolution_set_ = true;
    time_resolution_ = time_resolution;
  }

  void KMC_CoarseGrainSystem::initializeSystem(unordered_map<int, unordered_map<int, double>>& ratesOfAllSites) {

    LOG("Initializeing system", 1);

    if(!time_resolution_set_){
      throw runtime_error("You must first set the time resolution of the system "
          "before you can initialize the system.");
    }

    topology_->setRates(ratesOfAllSites);
		crude_.initializeSystem(ratesOfAllSites);
	}

  void KMC_CoarseGrainSystem::checkRates(unordered_map<int, unordered_map<int, double>>& ratesOfAllSites){
    vector<int> siteIds;
    for ( pair<int,unordered_map<int,double>> & site_neigh_rates : ratesOfAllSites){
      siteIds.push_back(site_neigh_rates.first);
      for ( pair<int,double> & neigh_rate : site_neigh_rates){
        siteIds.push_back(neigh_rate.first);
      } 
    }

    sort(siteIds.begin(),siteIds.end());
    siteIds.erase(unique(siteIds.begin(),siteIds.end()), siteIds.end());

    for ( size_t index = 0; index < siteIds.size();++index) {
      for( pair<int,double> & neigh_rate : ratesOfAllSites[siteIds.at(index)]){
        if(neigh_rate.first>siteIds.at(index)){
          if(!ratesOfAllSites.count(neigh_rate.first)){
            cerr << "No rate from " << neigh_rate.first << " to " << siteIds.at(index) << endl;
            throw runtime_error("ERROR detailed balance compromised not all rates are balanced.");  
          }
          if(!ratesOfAllSites.at(neigh_rate.first).count(siteIds.at(index))){
            cerr << "No rate from " << neigh_rate.first << " to " << siteIds.at(index) << endl;
            throw runtime_error("ERROR detailed balance compromised not all rates are balanced.");  
          } 
        }
      }
    }
  }

  int KMC_CoarseGrainSystem::getVisitFrequencyOfSite(int siteId){
    if(!rates_.count(siteId)){
      throw invalid_argument("Site is not stored in the coarse grained system you"
          " cannot retrieve it's visit frequency.");                            
    }                                                                           

    int visits = crude_.getVisitFrequencyOfSite(siteId);                
    if(topology_.partOfCluster(siteId)){                                           
      int cluster_id = topology_.getClusterIdOfSite(siteId);                       
      visits += topology_.getKMC_Cluster(cluster_id).getVisitFrequency(siteId); 
    }                                                                           
    return visits;  
    
  }

  void KMC_CoarseGrainSystem::initializeWalkers(vector<pair<int,KMC_Walker>>& walkers) {

    LOG("Initializeing walkers", 1);

		crude_.initializeWalkers(walkers);
  }

  void KMC_CoarseGrainSystem::setMinCoarseGrainIterationThreshold(int threshold_min) {
    LOG("Setting minimum coarse graining threshold", 1);
    iteration_threshold_min_ = threshold_min;
    iteration_threshold_ = threshold_min;
  }

  void KMC_CoarseGrainSystem::setRandomSeed(const unsigned long seed) {
    topology_->setRandomSeed(seed);
		crude_.setRandomSeed(seed);
  }

  void KMC_CoarseGrainSystem::removeWalkerFromSystem(pair<int,KMC_Walker>& walker) {
    removeWalkerFromSystem(walker.second);
  }

  void KMC_CoarseGrainSystem::removeWalkerFromSystem(KMC_Walker& walker) {
    LOG("Walker is being removed from system", 1);
		crude_.removeWalkerFromSystem(walker);
  }

  int KMC_CoarseGrainSystem::getClusterIdOfSite(int siteId) {
    return topology_->getClusterIdOfSite(siteId);
  }

	void runCrude(KMC_CoarseGrainSystem & CGSystem, int ,KMC_Walker & walker){
		CGSystem.crude_.hop(walker);	
	}

	void runCoarse(KMC_CoarseGrainSystem & CGSystem, int walker_id,KMC_Walker & walker){
		cout << "CGSystem hop" << endl;
		const int & siteId = walker.getIdOfSiteCurrentlyOccupying();
		const int & siteToHopToId = walker.getPotentialSite();
		KMC_TopologyFeature * feature = CGSystem.topology_->features[siteId].feature(*CGSystem.topology_,siteId);
		cout << "Potential site " << endl;

		if(!CGSystem.topology_->siteExist(siteToHopToId)){
			cout << "Creating site does not exist" << endl;
			CGSystem.topology_->features[siteToHopToId].feature(*CGSystem.topology_,siteToHopToId);
		}
		KMC_TopologyFeature * feature_to_hop_to = CGSystem.topology_->features[siteToHopToId].feature(*CGSystem.topology_,siteToHopToId);
		cout << "Potential site done" << endl;
    if(!feature_to_hop_to->isOccupied(siteToHopToId)){
			cout << "1" << endl;
      feature->vacate(siteId);
      feature_to_hop_to->occupy(siteToHopToId);

      walker.occupySite(siteToHopToId);
      walker.setDwellTime(feature_to_hop_to->getDwellTime(walker_id));
      walker.setPotentialSite(feature_to_hop_to->pickNewSiteId(walker_id));
    }else{
			cout << "2" << endl;
      feature->vacate(siteId);
      feature->occupy(siteId);

      walker.setDwellTime(feature->getDwellTime(walker_id));
      walker.setPotentialSite(feature->pickNewSiteId(walker_id));
    }

		cout << "CGSystem hop done" << endl;
	}

 void KMC_CoarseGrainSystem::hop(int walker_id,KMC_Walker& walker) {
    const int & siteId = walker.getIdOfSiteCurrentlyOccupying();
		const int & siteToHopToId = walker.getPotentialSite();
		site_funct_[siteId].run(*this,walker_id,walker);
    ++iteration_;
    if(iteration_ > iteration_threshold_){
      if(iteration_threshold_min_!=constants::inf_iterations){
        if(coarseGrain_(siteToHopToId)){
          iteration_threshold_ = iteration_threshold_min_;
        }else{
          iteration_threshold_*=2;
        }
				cout << "coarse grain done" << endl;
      }
      iteration_ = 0;
    }
  }

  /****************************************************************************
   * Internal Private Functions
   ****************************************************************************/

  bool KMC_CoarseGrainSystem::coarseGrain_(int siteId){

		cout << "Calling coarse grain method" << endl;
		vector<int> basin_site_ids;
		bool explore_basin = true;
		size_t basin_id_count = 0;
		size_t old_basin_id_count = 0;
		while(explore_basin){
			explore_basin = false;
			BasinExplorer basin_explorer;
			basin_site_ids = basin_explorer.findBasin(*topology_,siteId);
			basin_id_count = basin_site_ids.size();	
			if(basin_site_ids.size()>1 && basin_site_ids.size()<6){

				unordered_map<int,unordered_map<int,double>> external_rates = topology_->getExternallyConnectedRates(basin_site_ids);	
				unordered_map<int,unordered_map<int,double>>::iterator  iter = external_rates.begin();
				pair<int,double> fastest_rate{iter->second.begin()->first, (iter->second.begin()->second)};
				pair<int,double> second_fastest_rate = fastest_rate;

				for(;iter!=external_rates.end();++iter){
					for( pair<int,double> neigh_rates : iter->second){
						if(neigh_rates.second > fastest_rate.second){
							second_fastest_rate = fastest_rate;
							fastest_rate.first = neigh_rates.first;
							fastest_rate.second = neigh_rates.second;
						}else if(fastest_rate==second_fastest_rate){
							second_fastest_rate.first = neigh_rates.first;
							second_fastest_rate.second = neigh_rates.second;
						}
					}
				}

				// Determine if fastest rate is 3 orders of magnitude faster than second fastest rate
				if(fastest_rate.second/second_fastest_rate.second>1E4){
					// In turn check to see if the fastest rate off the site is back onto the cluster
					if(!topology_->siteExist(fastest_rate.first)){
						topology_->features[fastest_rate.first].feature(*topology_,fastest_rate.first);
					}
					unordered_map<int,double> rates_from_potential_neighbor = topology_->getRates(fastest_rate.first);		
					unordered_map<int,double>::iterator iter2 = rates_from_potential_neighbor.begin(); 
					pair<int,double> fastest_rate2(iter2->first,iter2->second);
					pair<int,double> second_fastest_rate2 = fastest_rate2;
					for(;iter2!=rates_from_potential_neighbor.end();++iter2){
						if(iter2->second > fastest_rate2.second){
							second_fastest_rate2 = fastest_rate2;
							fastest_rate2.first = iter2->first;
							fastest_rate2.second = iter2->second;
						}else if(fastest_rate2==second_fastest_rate2){
							second_fastest_rate2.first = iter2->first;
							second_fastest_rate2.second = iter2->second;
						}
					}

					if(fastest_rate2.second/second_fastest_rate2.second>1E4){
						if(find(basin_site_ids.begin(),basin_site_ids.end(),fastest_rate2.first)!=basin_site_ids.end()){
							//throw runtime_error("Consider expanding basin");
							if(old_basin_id_count!=basin_id_count){
								old_basin_id_count = basin_id_count;
								explore_basin = true;
								siteId = fastest_rate2.first;
								cout << "Explore true" << endl;
							}
						}
					}
				}
			}
		}
		cout << "Done with loop" << endl;
		if(basin_site_ids.size()>1){
			double internal_time_limit = getInternalTimeLimit_(basin_site_ids);

			cout << "Sites satisfy " << endl;
			if( sitesSatisfyEquilibriumCondition_(basin_site_ids, internal_time_limit) ){
				auto sites_and_clusters = getClustersOfSites(basin_site_ids);
				auto number_clusters = countUniqueClusters(sites_and_clusters);

			cout << "1" << endl;
				if(number_clusters==1 &&
						sites_and_clusters.begin()->second==constants::unassignedId)
				{
			cout << "2" << endl;
					createCluster_(basin_site_ids,internal_time_limit);
					return true;
				}else if(number_clusters!=1){
					// Joint clusters and sites to an existing cluster
			cout << "3" << endl;
					int favored_clusterId = getFavoredClusterId(sites_and_clusters);
					mergeSitesAndClusters_(sites_and_clusters,favored_clusterId);
					cout << "4" << endl;
					return true;
				}
			}
		}
		cout << "Done with coarse grain" << endl;
    return false;
  }

  size_t countUniqueClusters(const unordered_map<int,int> & sites_and_clusters){
    set<int> clusters;
    for(auto site_and_cluster : sites_and_clusters){
      clusters.insert(site_and_cluster.second);
    }
    return clusters.size();
  }

  int getFavoredClusterId(unordered_map<int,int> sites_and_clusters){
    int clusterId = constants::unassignedId;
    for(auto site_and_cluster : sites_and_clusters){
      if(site_and_cluster.second != constants::unassignedId){
        if(clusterId==constants::unassignedId || 
            site_and_cluster.second < clusterId){
          clusterId = site_and_cluster.second;
        }
      }
    }
    return clusterId;
  }

  // The first int is the site id the second int is the cluster id 
  unordered_map<int,int> KMC_CoarseGrainSystem::getClustersOfSites(const vector<int> & siteIds){
    return topology_->getClustersOfSites(siteIds); 
  }

  int KMC_CoarseGrainSystem::createCluster_(vector<int> siteIds, double internal_time_limit) {
    LOG("Creating cluster from vector of sites", 1);

    KMC_Cluster cluster;
    cluster.setConvergenceMethod(KMC_Cluster::Method::converge_by_tolerance);
    cluster.setConvergenceTolerance(0.001);
    
    cluster.addSites(siteIds);
    cluster.updateProbabilitiesAndTimeConstant();

    double cluster_time_const = cluster.getTimeConstant();
    // Cut the resolution in half from what it would otherwise be otherwise not worth doing
    double res = cluster_time_const/(2*internal_time_limit);
    double allowed_resolution = cluster_time_const/time_resolution_;
    double chosen_resolution = res;

    // The coarser the resolution is the better
    if(allowed_resolution <  chosen_resolution) chosen_resolution=allowed_resolution;

    if(chosen_resolution<2.0) chosen_resolution=2.0;
   
    cluster.setResolution(chosen_resolution);

    topology_->addKMC_Cluster(cluster);


    return cluster.getId();
  }

  void KMC_CoarseGrainSystem::mergeSitesAndClusters_( unordered_map<int,int> sites_and_clusters,int favoredClusterId) {

    LOG("Merging sites to cluster", 1);
    vector<int> isolated_sites;
    unordered_set<int> cluster_ids;

    for (const pair<int,int> & site_and_cluster : sites_and_clusters) { 
				cout << "T0" << endl;
      if(site_and_cluster.second != favoredClusterId){ 
				cout << "T01" << endl;
        if (site_and_cluster.second == constants::unassignedId) {
          isolated_sites.push_back(site_and_cluster.first);
        } else {
          cluster_ids.insert(site_and_cluster.second);
        }
        topology_->features[site_and_cluster.first].feature = returnCluster; 
				cout << "T1" << endl;
				topology_->setSitesClusterId(site_and_cluster.first,favoredClusterId);
				cout << "T2" << endl;
      }

				cout << "T3" << endl;
			site_funct_[site_and_cluster.first].run = runCoarse;
				cout << "T4" << endl;
    }
				cout << "T5" << endl;
    topology_->getKMC_Cluster(favoredClusterId).addSites(isolated_sites);
    topology_->getKMC_Cluster(favoredClusterId).updateProbabilitiesAndTimeConstant();
				cout << "T6" << endl;
    for(auto clusterId : cluster_ids ){
				cout << "T7" << endl;
      topology_->mergeClusters(favoredClusterId,clusterId);
				cout << "T8" << endl;
    }

  }

  double KMC_CoarseGrainSystem::getExternalTimeLimit_(const vector<int> & siteIds ){
    LOG("Getting the external time limit of a cluster", 1);
    unordered_set<int> internal_sites;
    for(const int & site_id : siteIds){
      internal_sites.insert(site_id);
    }

    double max_rate_off = 0; 
    for(const int & site_id : siteIds){
      const unordered_map<int,double> & neigh_and_rates = rates_->at(site_id);
      for( const pair<int,double> & neigh_and_rate : neigh_and_rates){
        if(internal_sites.count(neigh_and_rate.first)==0){
          if((neigh_and_rate.second) > max_rate_off){
            max_rate_off = (neigh_and_rate.second);
          }
        }
      }
    }
    return max_rate_off;
  }

double KMC_CoarseGrainSystem::getInternalTimeLimit_(vector<int> siteIds ){
  LOG("Getting the internal time limit of a cluster", 1);

  auto nodes = convertSitesToEmptySharedNodes(siteIds);
  unordered_map<int, weak_ptr<GraphNode<string>>> nodes_weak;
  for (auto node_iter : nodes) nodes_weak[node_iter.first] = node_iter.second;
  auto edges = convertSitesOutgoingRatesToTimeSharedWeightedEdges<vector<shared_ptr<Edge>>>(
      *topology_,
      siteIds);

  list<weak_ptr<Edge>> edges_weak(edges.begin(), edges.end());

  auto graph_ptr = Graph<string>(edges_weak, nodes_weak);

  unordered_map<pair<int, int>, double,hash_functions::hash> verticesAndtimes =
      maxMinimumDistanceBetweenEveryVertex<string>(graph_ptr);

  double maxtime = 0.0;
  for (auto verticesAndTime : verticesAndtimes) {
    if (verticesAndTime.second > maxtime) maxtime = verticesAndTime.second;
  }
  return maxtime;
}

// Its not worth creating a cluster unless the time is at least cut in half
// And it is not allowed if the sample time is smaller than than the simulated
// time of the cluster. The cluster has to be updated at a minimum once between
// each measurment (time_resolution). If this is not done the noise will not
// correctly show up in the data.  
// The number 25 is the ratio needed between hops within the cluster to hops
// outside of the cluster in order to see performance gains.
bool KMC_CoarseGrainSystem::sitesSatisfyEquilibriumCondition_(
    vector<int> siteIds, double maxtime) {

  LOG("Checking if sites satisfy equilibrium condition", 1);
  double timeConstant = getTimeConstantFromSitesToNeighbors_(siteIds);
  double time_to_traverse_cluster = maxtime*minimum_coarse_graining_resolution_;
  return timeConstant > time_to_traverse_cluster*performance_ratio_ && time_to_traverse_cluster< time_resolution_;// && ratio>25;
}

double KMC_CoarseGrainSystem::getTimeConstantFromSitesToNeighbors_(
   const vector<int> & siteIds) const {

  return topology_->getTimeConstantFromSitesToNeighbors(siteIds);
}

unordered_map<int,vector<int>> KMC_CoarseGrainSystem::getClusters(){
  return topology_->getClusters();
}

unordered_map<int,double> KMC_CoarseGrainSystem::getResolutionOfClusters(){
  return topology_->getResolutionOfClusters();
}

unordered_map<int,double> KMC_CoarseGrainSystem::getTimeIncrementOfClusters(){
  return topology_->getTimeIncrementOfClusters();
}

int KMC_CoarseGrainSystem::getFavoredClusterId_(vector<int> siteIds) {
  return topology_->getFavoredClusterId(siteIds);
}


}
