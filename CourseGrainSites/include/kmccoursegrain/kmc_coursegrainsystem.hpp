#ifndef KMCCOURSEGRAIN_KMC_COURSEGRAINSYSTEM_HPP
#define KMCCOURSEGRAIN_KMC_COURSEGRAINSYSTEM_HPP

#include <map>
#include <unordered_set>
#include <unordered_map>
#include <memory>
#include <vector>

#include "../../src/libkmccoursegrain/kmc_site_container.hpp"
#include "../../src/libkmccoursegrain/topologyfeatures/kmc_site.hpp"
#include "../../src/libkmccoursegrain/topologyfeatures/kmc_cluster.hpp"
#include "../../src/libkmccoursegrain/topologyfeatures/kmc_topology_feature.hpp"

#include "kmc_constants.hpp"

namespace ugly {
template <typename... Ts>
class Graph;
}

namespace kmccoursegrain {

class KMC_Particle;

/**
 * \brief Course Grain System allows abstraction of renormalization of sites
 *
 * This class simply uses pointers to the relevant hop rates. It will simulate a
 * particle hopping through the system of sites. If a large number of compute
 * cycles are expended moving a particle between two low energy sites they will
 * be course grained, or renormalized so that the probabilities and time spent
 * on each site will be the same but the number of compute cycles will be
 * significantly reduced.
 **/
class KMC_CourseGrainSystem {

 public:
  /**
   * \brief Constuctor for course grained system
   *
   * The constructor by default seeds the random number generator based on the
   * current time. Furthermore, course graining of two sites is by default set
   * to a threshold of 20. Meaning a particle must remember moving back and
   * forth between two sites at least 20 times before the sites are course
   * grained.
   **/
  KMC_CourseGrainSystem()
      : seed_set_(false),
        max_cluster_resolution_(20),
        minimum_course_graining_resolution_(2),
        iteration_threshold_(10000),
        iteration_threshold_min_(10000){};

  /**
   * \brief This will correctly initialize the system
   *
   * Essentially pointers to all the rates are stored in this class. Pointers
   * are used in case any of the rates are changed, there will be no need to
   * pass them back into the couse grain object. An updated function could
   * simply be called. This function must be called before the particles can
   * be initialized `initializeParticles` and before a hopping event is called
   * on a particle `hop`.
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
  void initializeSystem(std::unordered_map<int, std::unordered_map<int, double*>> ratesOfAllSites);

  /**
   * \brief Initialize particle dwell times and future hop site id
   *
   * This function can only be called after initializing the sites. Each
   * particle mst have also been placed on a site in the system. As in if
   * there are sites 1 2 and 3. Then the particles must exist on at least one
   * of these sites before they are passed in. The function will then update
   * their dwell times as well as providing a potential future hopping site.
   *
   * \param[in] particles a vector of pointers to the particles
   **/
  void initializeParticles(std::vector<KMC_Particle>& particles);

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
   * \brief Make the particle hop to a site in the system
   *
   * Once a particle has been initialized, i.e. it has a dwell time it is
   * located on a site in the system and it has a stored potential site it
   * will be hopping to. This function will be called on it. It will move the
   * particle if necessary and will course grain sites/renormalize sites if
   * necessary.
   *
   * \param[in] particle
   **/
  void hop(KMC_Particle& particle);

  /**
   * \brief Remove the particle from the system
   **/
  void removeParticleFromSystem(KMC_Particle& particle);

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
   * \brief Determines how often to check for course graining
   *
   * This function determines the number of iterations that will pass before 
   * the algorithm will check for course graining 
   *
   * \param[in] thershold
   **/
  void setMinCourseGrainIterationThreshold(int threshold_min);
  int getMinCourseGrainIterationThreshold();

  std::vector<std::vector<int>> getClusters();

  /**
   * \brief set and get the course graining resolution
   *
   * The smaller this value is the faster the algorithm should be. However, this
   * is at the cost of reproducing the noise of the simulation. If you are 
   * wanting to capture the noise this needs to be higher. However, it is lower 
   * the more effiecient the course graining should be.
   **/
  int getMaxCourseGrainResolution() { return max_cluster_resolution_; }
  void setMaxCourseGrainResolution(int max_cluster_resolution) {
    max_cluster_resolution_ = max_cluster_resolution;
  }

 private:
  /// Depicts whether a random seed has been set, to yield reproducable data
  bool seed_set_;

  /// The random seed
  unsigned long seed_;

  /// The resolution of the clusters. Essentially how many hops will a particle
  /// move within the cluster before it is likely to leave, the point of this
  /// is to at least to a small degree conserve the noise.
  int max_cluster_resolution_;

  /// This should be set to a value of 2, it is used to determine if course 
  /// graining should occur. If the time to hop off the potential sites in
  /// the course grained cluster is less than twice as long it is not worth
  /// course graining. 
  int minimum_course_graining_resolution_;

  /// Keeps track of the number of iterations. Is reset after passing the
  /// iteration threshold. 
  int iteration_;

  /// How many interactions occur before course graining is tested
  int iteration_threshold_;

  /// The iteration threshold is reset to the min value if a cluster is found
  int iteration_threshold_min_;

  std::unordered_map<int, KMC_TopologyFeature *> topology_features_;
  /// Stores smart pointers to all the sites
  KMC_Site_Container sites_;

  /// Stores smart pointers to all the clusters
  std::unordered_map<int, KMC_Cluster> clusters_;

  void courseGrainSiteIfNeeded_(KMC_Particle& particle);

  /**
   * \brief Determines if it is appropriate to coursegrain the sites
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

  bool courseGrain_(int siteId);
  std::unordered_map<int,int> getClustersOfSites(std::vector<int> siteIds);
  int createCluster_(std::vector<int> siteIds,double internal_time_limit);
  void mergeSitesAndClusters_(std::unordered_map<int,int> sites_and_clusters, int clusterId);
  double getMinimumTimeConstantFromSitesToNeighbors_(std::vector<int> siteIds);
  std::unordered_map<int,double> filterSites_();
  std::vector<std::vector<int>> breakIntoIslands_(std::unordered_map<int,double> relevant_sites);
};
}
#endif  // KMCCOURSEGRAIN_KMC_COURSEGRAINSYSTEM_HPP
