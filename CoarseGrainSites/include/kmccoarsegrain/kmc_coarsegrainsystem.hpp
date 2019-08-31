#ifndef KMCCOARSEGRAIN_KMC_COARSEGRAINSYSTEM_HPP
#define KMCCOARSEGRAIN_KMC_COARSEGRAINSYSTEM_HPP

#include <map>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <vector>

#include "kmc_constants.hpp"
#include "kmc_crude.hpp"

namespace ugly {
template <typename... Ts>
class Graph;
}

namespace kmccoarsegrain {

class KMC_CoarseGrainSystem; 
class KMC_Dynamic_Topology;
class KMC_TopologyFeature;
class KMC_Walker;


	void runCrude(KMC_CoarseGrainSystem & CGsystem,int walker_id,KMC_Walker & walker);
	void runCoarse(KMC_CoarseGrainSystem & CGsystem,int walker_id,KMC_Walker & walker);
/**
 * \brief Coarse Grain System allows abstraction of renormalization of sites
 *
 * This class simply uses pointers to the relevant hop rates. It will simulate a
 * walker hopping through the system of sites. If a large number of compute
 * cycles are expended moving a walker between two low energy sites they will
 * be coarse grained, or renormalized so that the probabilities and time spent
 * on each site will be the same but the number of compute cycles will be
 * significantly reduced.
 **/
class KMC_CoarseGrainSystem {

 public:
  /**
   * \brief Constuctor for coarse grained system
   *
   * The constructor by default seeds the random number generator based on the
   * current time. Furthermore, coarse graining of two sites is by default set
   * to a threshold of 20. Meaning a walker must remember moving back and
   * forth between two sites at least 20 times before the sites are coarse
   * grained.
   **/
  KMC_CoarseGrainSystem();
  ~KMC_CoarseGrainSystem();
  /**
   * \brief This will correctly initialize the system
   *
   * Essentially pointers to all the rates are stored in this class. Pointers
   * are used in case any of the rates are changed, there will be no need to
   * pass them back into the couse grain object. An updated function could
   * simply be called. This function must be called before the walkers can
   * be initialized `initializeWalkerss` and before a hopping event is called
   * on a walker `hop`.
   *
   * \param[in] ratesOfAllSites this is a map of maps the first int is site i
   * the value of which is a second map of sites j which are all neighbors of
   * site i. The double of the final map is a pointer to the actual rate.
   * Consider the following to understand how the maps are structured. I have
   * a site 1 which has three neighbors and a site 2 which has 2 neighbors
   *
   * site4 - site1 - site2 - site3
   *           |
   *         site5
   *
   * Each line - is composed of two rates <- and ->. So if I were to store the
   * rates in the map it would look like this
   *
   * map<int,map<int,double *>> rates;
   *
   * // Rates from site 1
   *
   * rates[1][2] = &rateFrom1to2;
   * rates[1][4] = &rateFrom1to4;
   * rates[1][5] = &rateFrom1to5;
   *
   * // Rates from site 2
   *
   * rates[2][1] = &rateFrom2to1;
   * rates[2][3] = &rateFrom2to3;
   *
   * // Rates from site 3
   *
   * rates[3][2] = &rateFrom3to2;
   *
   * // Rates from site 4
   *
   * rates[4][1] = &rateFrom4to1;
   *
   * // Rates from site 5
   *
   * rates[5][1] = &rateFrom5to1;
   *
   * Where each of the rateFrom variables is a double
   **/
  void initializeSystem(std::unordered_map<int, std::unordered_map<int, double>> &ratesOfAllSites);

  const std::unordered_map<int,std::unordered_map<int,double>> rates();

  /**
   * @brief Checks that detailed balance is satisfied
   *
   * Every rate off a site to a neighbor must be balanced by a rate from
   * the neighbor to the site
   *
   * @param ratesOfAllSites
   */
  void checkRates(std::unordered_map<int, std::unordered_map<int, double>> &ratesOfAllSites);

  /**
   * \brief Initialize walker dwell times and future hop site id
   *
   * This function can only be called after initializing the sites. Each
   * walker mst have also been placed on a site in the system. As in if
   * there are sites 1 2 and 3. Then the walkers must exist on at least one
   * of these sites before they are passed in. The function will then update
   * their dwell times as well as providing a potential future hopping site.
   *
   * \param[in] walkers a vector of pointers to the walkers
   **/
  void initializeWalkers(std::vector<std::pair<int,KMC_Walker>>& walkers);

  /**
   * \brief Define the seed for the random number generator
   *
   * This allows the user to create reproducable results if desired. By
   * default the seed will be determined from the time.
   *
   * \param[in] seed
   **/
  void setRandomSeed(const unsigned long seed);

  /**
   * \brief Make the walker hop to a site in the system
   *
   * Once a walker has been initialized, i.e. it has a dwell time it is
   * located on a site in the system and it has a stored potential site it
   * will be hopping to. This function will be called on it. It will move the
   * walker if necessary and will coarse grain sites/renormalize sites if
   * necessary.
   *
   * \param[in] walker
   **/
  void hop(std::pair<const int, KMC_Walker>& walker);
  void hop(int walker_id, KMC_Walker& walker);
  //void hop(KMC_Walker& walker);

  /**
   * \brief Remove the walker from the system
   **/
  void removeWalkerFromSystem(std::pair<int,KMC_Walker>& walker);
  void removeWalkerFromSystem(KMC_Walker& walker);

  /**
   * \brief Determine if the site is part of a cluster
   *
   * Will return the id of the cluster the site is a part of or if the site
   * is not part of the cluster it will return constants::unassignedId
   *
   * \param[in] siteId
   *
   * \return clusterId
   **/
  int getClusterIdOfSite(int siteId);

  int getVisitFrequencyOfSite(int siteId);

  /**
   * \brief Determines how often to check for coarse graining
   *
   * This function determines the number of iterations that will pass before 
   * the algorithm will check for coarse graining. If one desires to turn the
   * cluster algorithm off the threshold can be set with a value of:
   * constants::inf_iterations. 
   *
   * \param[in] thershold
   **/
  void setMinCoarseGrainIterationThreshold(int threshold_min);
  int getMinCoarseGrainIterationThreshold();


  /**
   * @brief Return the clusters
   *
   * vector of the indices of the cluster beginning with the id of the cluster
   *
   * @return 
   */
  std::unordered_map<int,std::vector<int>> getClusters();

  std::unordered_map<int,double> getResolutionOfClusters();
  std::unordered_map<int,double> getTimeIncrementOfClusters();
  /**
   * \brief Determines how fine grained the time is allowed to be
   *
   * When running simulations where you are only sampling the simulation every
   * few time increments, a fine grained resolution may not be needed. Reducing
   * the resolution allows the cluster algorithm to perform better.
   **/
  double getTimeResolution();
  void setTimeResolution(double time_resolution);

  /**
   * @brief Adjusts how easy it is to create a cluster.
   *
   * This parameter should not be changed unless you know what you are doing. 
   * Essentially, the larger it is, delays when the cluster will satisfy the
   * coarse graining criteria. 
   *
   * When determining if sites will make a good cluster, the rates off the 
   * potential cluster sites are used to create a time constant. If this time
   * constant is greater than the time it takes for a particle to move within
   * the cluster it would normally satify the mathematical arguments for 
   * coarse graining. However, there is extra overhead associated with using
   * clusters, the overhead becomes negligable the larger the different between
   * the time constant and the internal time. 
   *
   * The parameter essentially delays when coarse graining will occur so that
   * the algorithm does not take a performance penalty be coarse graining sites
   * that do not have a big enough different between the time constant and the
   * internal time of crossing. 
   *
   * @param double
   */
  void setPerformanceRatio(double performance_ratio) { 
    performance_ratio_ = performance_ratio;
  }
 private:

	KMC_Crude crude_;

	// Crude implementation only stores the cummulitive probability dist but
	// we need the actual probs in order coarse grain the clusters
	std::unordered_map<int,std::unordered_map<int,double>> site_neigh_prob;
  /// Performance ratio
  double performance_ratio_;
/*
//  std::unordered_map<int,std::unordered_map<int,double>> * rates_;
  /// Depicts whether a random seed has been set, to yield reproducable data
  bool seed_set_;

  /// The random seed
  unsigned long seed_;
*/
  bool time_resolution_set_;
  /// The resolution of the clusters. Essentially how many hops will a walker
  /// move within the cluster before it is likely to leave, the point of this
  /// is to at least to a small degree conserve the noise.
  //int max_cluster_resolution_;
  double time_resolution_;

  /// This should be set to a value of 2, it is used to determine if coarse 
  /// graining should occur. If the time to hop off the potential sites in
  /// the coarse grained cluster is less than twice as long it is not worth
  /// coarse graining. 
  int minimum_coarse_graining_resolution_;

  /// Keeps track of the number of iterations. Is reset after passing the
  /// iteration threshold. 
  int iteration_;

  /// How many interactions occur before coarse graining is tested
  int iteration_threshold_;

  /// The iteration threshold is reset to the min value if a cluster is found
  int iteration_threshold_min_;

  /// Must be a pointer so that we do not have to include private header
  std::unique_ptr<KMC_Dynamic_Topology> topology_;

	void coarseGrainSiteIfNeeded_(KMC_Walker& walker);

  /**
   * \brief Determines if it is appropriate to coarsegrain the sites
   *
   * This function looks to see if the Markov property holds for the sites of
   * interest. This is done by creating a graph consisting of the sites that
   * will be placed in the cluster. Once this is done we determine the fastest
   * time it takes to cross from one side of the cluster to the other. If this
   * time is much less than the transition of the sites in the cluster we know
   * that we can approximate it as being in a local equilibrium
   *
   * \param[in] siteIds - vector of site ids that will make up the cluster
   *
   * \return true if the sites satisfy the condition false otherwise
   **/
  bool sitesSatisfyEquilibriumCondition_(std::vector<int> siteIds, double maxtime);

  double getInternalTimeLimit_(std::vector<int> siteIds);

  /**
   * @brief Gets the fastest rate off the basin sites
   *
   * @param siteIds
   *
   * @return 
   */
  double getExternalTimeLimit_(const std::vector<int> & siteIds);

  /**
   * \brief Determines that the cluster id should be if sites are to be merged
   *
   * If there are several sites that could make a cluster we determine if they
   * are already part of a cluster or not. If they are part of a cluster or more
   * than one cluster is found. We will merge clusters to the cluster will the
   * smallest cluster id.
   *
   * \param[in] siteIds site ids that will potentially make up a cluster
   *
   * \return int value that represents the smallest clsuter id or else it
   * returns constant::unassignedId
   **/
  int getFavoredClusterId_(std::vector<int> siteIds);


	struct DefaultSiteFunction {
		void (*run)(KMC_CoarseGrainSystem & CGSystem, int walker_id,KMC_Walker & walker) = runCrude;
	};

	friend class KMC_Cluster;
	friend void runCrude(KMC_CoarseGrainSystem & CGsystem,int walker_id,KMC_Walker & walker);
	friend void runCoarse(KMC_CoarseGrainSystem & CGsystem,int walker_id,KMC_Walker & walker);
	std::unordered_map<int,DefaultSiteFunction> site_funct_;


  bool coarseGrain_(int siteId);
  std::unordered_map<int,int> getClustersOfSites(const std::vector<int> & siteIds);
  int createCluster_(std::vector<int> siteIds,double internal_time_limit);
  void mergeSitesAndClusters_(std::unordered_map<int,int> sites_and_clusters, int clusterId);
  double getTimeConstantFromSitesToNeighbors_(const std::vector<int> & siteIds) const;
  std::unordered_map<int,double> filterSites_();
  std::vector<std::vector<int>> breakIntoIslands_(std::unordered_map<int,double> relevant_sites);
};
}
#endif  // KMCCOARSEGRAIN_KMC_COARSEGRAINSYSTEM_HPP
