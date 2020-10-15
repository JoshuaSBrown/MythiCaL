#ifndef MYTHICAL_CLUSTER_HPP
#define MYTHICAL_CLUSTER_HPP

#include <cassert>
#include <list>
#include <memory>
#include <random>
#include <vector>
#include <unordered_map>

#include "topology_feature.hpp"
#include "site.hpp"

namespace mythical {

/**
 * \brief Coarse graining of sites is handled by the Cluster class
 *
 * This class is responsible for coarse graining problem sites that occur
 * during kinetic monte carlo simulations. It works by assuming that the sites
 * contained in the cluster are in equilibrium. The master equation is then
 * solved in an interative procedure until the probabilities of a charge
 * exisiting on each of the sites can be determined. Once this is known the
 * probability of hopping to sites surrounding the cluster can also be
 * calculated. As well as the dwell time. etc...
 **/
class Cluster : public TopologyFeature {

  public:
  /**
   * \brief Constructor for cluster
   *
   * Specific default values are associated with the cluster including the
   * following:
   *
   * iterations: 3 - Number of iterations to reach convergence
   *
   * convergence method: iterations*(number of sites)
   *
   * resolution: 20 - The maximum amount of time simulated before a hopping
   * event
   *
   * random number generator seed: based on time
   *
   * It is also worth noting that each cluster will be given a different id
   * on creation. This is automated in the constructor. Such that every new
   * cluster will be given an integer id greater than the last one created.
   **/
  Cluster();

  /**
   * \brief Convergence Methods
   *
   * The convergence methods define how the convergence is handled, for
   * sites that are considered to be part of a cluster. There are two
   * parameters that are associated with these methods the number of
   * iterations and the tolerance.
   *
   * converge_by_iterations_per_cluster
   *
   * This method relies only on the number of iterations that have been set.
   * it works by using a fixed set of iterations to establish convergence of
   * the master equation approach used in the convergence technique. So if
   * this method is chosen and you have also set the number of iterations to
   * 10. Only 10 iterations will be used regardless of the number of sites
   * in the cluster.
   *
   * converge_by_iterations_per_site
   *
   * This method only requires that the number of iterations be specified.
   * The number of iterations used for the convergence of the Master
   * Equation of the sites in the cluster is calculated by multiplying the
   * number of iterations by the number of sites in the cluster. Thus if 3
   * iterations are specified and there are 5 sites in the cluster a total
   * of 15 iterations will be used.
   *
   * converge_by_tolerance
   *
   * This method only requires that the tolerance is specified. If it is
   * chosen convergence of the master equation continues for an unspecified
   * number of iterations but until the maximum difference of the sites
   * probabilities is less than the tolerance.
   **/
  enum Method {
    converge_by_iterations_per_cluster,
    converge_by_iterations_per_site,
    converge_by_tolerance
  };

  /**
   * \brief Adds a site to a cluster
   *
   * This function adds a site/sites to a cluster. An arbitrary number of
   * sites may be added. However, an error will be thrown if you attempt to
   * add the same site more than once.
   *
   * \param[in] site a shared pointer to a site
   **/
  void addSite(Site& site);
  void addSites(std::vector<Site>& sites);

  /**
   * \brief will update the probabilities and time constant stored in the
   * cluster
   *
   * Note that this function is needed if a site stored in the cluster is
   * changed. The cluster will have no idea if the site is changed unless
   * the update function is called. So for instance if site1 had a rate of
   * 500 to its neighbor when it is added to the cluster but later this rate
   * is changed to 0.04 the update must be called.
   **/
  void updateProbabilitiesAndTimeConstant();

  /**
   * \brief Determines if a site is the cluster
   *
   * \param[in] siteId integer value representing the id of the site
   *
   * \return True if the site is in the cluster false otherwise
   **/
  bool siteIsInCluster(const int siteId) const {
    return sitesInCluster_.count(siteId);
  }

  /**
   * \brief Create a vector storing all the ids of the sites that are in the
   * Cluster
   *
   * \return A vector of shared pointers to the sites
   **/
  std::vector<Site> getSitesInCluster() const;

  std::vector<int> getSiteIdsInCluster() const;
  std::vector<int> getSiteIdsNeighboringCluster() const;

  /**
   * \brief Returns the number of sites in the cluster
   **/
  int getNumberOfSitesInCluster() const { return sitesInCluster_.size(); }

  /**
   * \brief Calculates the probability of hopping to a site in the cluster
   *
   * If the site id passed in does not correspond to a site that is in the
   * cluster this function will throw an error.
   *
   * \param[in] siteId integer value indicating the sites id
   *
   * \return a double corresponding to the probability
   **/
  double getProbabilityOfOccupyingInternalSite(const int siteId);

  /**
   * \brief Move the sites in one cluster to another
   *
   * Here we are moving the sites from one cluster to the other. After doing
   * this any internal values that need to be recalibrated are.
   *
   * \param[in] cluster a smart pointer to a cluster
   **/
  void migrateSitesFrom(Cluster& cluster);

  /**
   * \brief Set the resolution of the cluster
   *
   * The resolution is an important parameter defining how large of a time
   * step the cluster will use. The larger the resolution the smaller the
   * improvement in computational performance. However, the smaller it is
   * the more coarse the represenation of the approximation. I have found
   * that values between 20-50 are a good starting point. The maximum dwell
   * time is defined as:
   *
   * max(dwell) = timeConstant/resolution
   *
   * \param[in] resolution integer value defining how coarse grained the
   * cluster is
   **/
  void setResolution(const double resolution) { 
    // Must be greater than 1
    assert(resolution>=2.0); 
    resolution_ = resolution; 
    time_increment_ = TopologyFeature::escape_time_constant_/static_cast<double>(resolution_);
  }

  double getResolution() const {
    return resolution_;
  }

  double getTimeIncrement() const {
    return time_increment_;
  }
  /**
   * \brief Pick the next site a particle will hop too
   *
   * This is one of the core methods. Will essentially pick a site within or
   * neighboring the cluster based on the calculated probabilities and
   * return the site id.
   *
   * \return site id generated to reproduce the probability of a particle
   * moving to it
   **/
  int pickNewSiteId(const int & walker_id);
  //int pickNewSiteId();

  /**
   * \brief Set the convergence method
   *
   * A detailed description of the methods is avaible above and is an enum.
   *
   * \param[in] convergence_method this is an enum determining how
   * convergence of the Master Equation should be handled
   **/
  void setConvergenceMethod(const Method convergence_method) {
    convergence_method_ = convergence_method;
  }

  /**
   * \brief Sets the convergence tolerance of the Master Equation approach
   *
   * This is only needed if the method specified is converge_by_tolerance. A
   * default value of 0.01 will have already been specified.
   *
   * \param[in] tolerance
   **/
  void setConvergenceTolerance(double tolerance);

  /**
   * \brief Return the convergence tolerance
   **/
  double getConvergenceTolerance() const { return convergenceTolerance_; }

  /**
   * \brief Sets the number of iterations used for converging the master
   * Equation
   *
   * The methods converge_by_iterations_per_cluster and
   * converge_by_iterations_per_site are the only ones that use this
   * parameter.
   **/
  void setConvergenceIterations(const long iterations);

  /**
   * \brief Return the number of convergence iterations
   **/
  long getConvergenceIterations() const { return iterations_; }

  /**
   * \brief Returns a probability of a particle moving to a neighbor
   *
   * \param[in] neighId the site id of the neighbor
   *
   * \return a probability
   **/
  double getProbabilityOfHoppingToNeighborOfCluster(const int neighId);

  /**
   * \brief Returns the dwell time, each call will return a different value
   **/
  double getDwellTime(const int & walker_id);
  //double getDwellTime();

  double getFastestRateOffCluster();

  void setVisitFrequency(int frequency,const int & siteId);
  int getVisitFrequency(const int & siteId);

  /**
   * \brief Prints the contents of the cluster
   **/
  friend std::ostream& operator<<(std::ostream& os,
                                  const Cluster& cluster);

 private:
  /************************************************************************
   * Local Cluster Variables
   ************************************************************************/
  /// Previous number of cluster visits
  int prev_total_visit_freq_;

  /// Relates to how coarse grained the dwell time will be
  double resolution_;

  /// Number of iterations used to solve the master equation
  long iterations_;

  /// Tolerance used to determine when the master equation has been solved
  double convergenceTolerance_;

  /// Type of convergence used to solve the master equation
  Method convergence_method_;

  /// Time increment of the cluster
  double time_increment_;

  double internal_time_constant_;
  /**
   * \brief Stores the particle and its total dwell time on the cluster
   *
   * This is an internal map, the dwell time returned by getDwellTime will 
   * return the actual dwell time if the resolution were 1. Whenever the charge
   * has been on the cluster longer than the dwell time stored by the cluster
   * it will then be able to hop to a site external to the cluster. 
   **/
  std::unordered_map<int,double> remaining_walker_dwell_times_;

  /**
   * \brief Stores the probability of hopping to each of the neighbors
   *
   * The first integer in the pair is the neighbor site id, the double is a
   * probability given as a value between 0 and 1.
   **/
  std::vector<std::pair<int, double>> probabilityHopToNeighbor_;
  std::vector<std::pair<int, double>> cumulitive_probabilityHopToNeighbor_;

  /**
   * \brief Stores the internal dwell time of the sites in the cluster
   *
   * In other words this stores the dwell times as if there are no rates to
   * sites neighboring the cluster. 
   **/
  std::unordered_map<int,double> internal_dwell_time_;

  /**
   * \brief Stores the number of times a site in the cluster is visited
   *
   **/
  std::unordered_map<int,double> site_visits_;

  /**
   * \brief Stores the rates from each site in the cluster to the neighbors
   *
   * The first int is the site id of an internal site, the double is a sum of
   * all the rates going from the internal site to sites neighboring the cluster.
   **/
  std::unordered_map<int, double> sumOfEscapeRateFromSiteToNeighbor_;
  /**
   * Same as above but the double is the sum of the rates going to other sites
   * within the cluster. 
   **/
  std::unordered_map<int, double> sumOfEscapeRateFromSiteToInternalSite_;

  /**
   * \brief Stores the pointers to sites that are in the cluster
   **/
  std::unordered_map<int, Site> sitesInCluster_;

  std::unordered_map<int,double> probabilityHopOffInternalSite_;
  std::unordered_map<int,double> probabilityHopBetweenInternalSite_;

  /**
   * \brief The probability of a particle being on each of the sites
   *
   * The first int is the site id of a site within the cluster, the double
   * is a probability between 0 and 1 which is found from solving the
   * Master Eqaution.
   **/
  std::unordered_map<int, double> probabilityOnSite_;

  std::vector<std::pair<int,double>> probabilityHopToInternalSite_;
  std::vector<std::pair<int,double>> cumulitive_probabilityHopToInternalSite_;

  /************************************************************************
   * Local Cluster Functions
   ************************************************************************/

  //void incrementVisits_();

  std::unordered_map<int,int> getVisitFrequencies_();

  /// Will solve the Master Equation
  void solveMasterEquation_();

  /**
   * \brief Picks a neighboring site of the cluster
   *
   * If it is determined that the particle will be moving out of the cluster
   * this function will with the correct probability distribution choose one
   * of the neighbors.
   *
   * \return the site id of one of the neighbors
   **/
  int pickClusterNeighbor_(const int & walker_id);
  //int pickClusterNeighbor_();

  /**
   * \brief Picks a site within the cluster
   *
   * If it is determined that a particle will not escape the cluster this
   * function using the correct probability distribution will pick a site
   * within the cluster to move to.
   *
   * \return the site id of a site within the cluster
   **/
  int pickInternalSite_();

  /**
     * \brief Calculates the time constant used to calculate the dwell time
     *
     * This time constant determines how much time it will take for the
     * particle to escape the cluster. It is essentially a constant in an
     * exponent describing the escape rate. E.g.
     *
     * f(t) = A exp(-beta * t)
     *
     * A is some constant, t is time, and f(x) describes the probability of
     * escaping the cluster and where:
     *
     * beta = 1/tau
     *
     * And tau is the time constant calculated here. So the larger tau is the
     * longer it will take for the particle to escape.
     **/
    void calculateEscapeTimeConstant_();
    void calculateInternalTimeConstant_();
    /**
     * \brief Determines if a particle will stay within a cluster or not
     *
     * This function uses a random number to determine if a particle will stay
     * within a cluster or hop to a site neighboring it.
     *
     * \return True if the particle should stay within the cluster, False if
     * it should not
     **/
    bool hopWithinCluster_(const int & walker_id) const;
    //bool hopWithinCluster_();

    // First int is the Id of a site within the cluster
    // pair - first int is the id of the site neighboring the cluster
    // double is the rate
    std::unordered_map<int, std::unordered_map<int, double>>
        getRatesToNeighborsOfCluster_();

    std::unordered_map<int, std::unordered_map<int, double>>
        getRatesBetweenInternalSites_();

    void iterate_();

    void calculateProbabilityHopToNeighbors_();
    void calculateProbabilityHopToInternalSite_();
    void calculateProbabilityHopOffInternalSite_();
    void calculateProbabilityHopBetweenInternalSite_();
    void calculateInternalDwellTimes_();
    void calculateSumOfEscapeRatesFromSitesToTheirNeighbors_();
    void calculateSumOfEscapeRatesFromSitesToInternalSites_();

    void initializeProbabilityOnSites_();

    /**
     * \brief Will grab all the internal rates going to each site in the
     * cluster
     *
     * This function organises the rates describing movement within a cluster
     * such that the first integer in the map is the internal neighbor id of
     * a site. The vector of pairs contains the directional rate going from
     * the neighbor to a site. Thus the integer in the pair is the site id
     * and the double is the rate.
     *
     * E.g.
     *
     * site1 <- internal_neigh1 -> site2
     *                |
     *                v
     *              site3
     *
     * The diagram explains how the contents are stored in the map.
     **/
    std::unordered_map<int, std::vector<std::pair<int, double>>>
        getInternalRatesFromNeighborsComingToSite_();

    friend void occupyCluster_(TopologyFeature*,const int&);
    friend void vacateCluster_(TopologyFeature*,const int&);
    friend bool isOccupiedCluster_(TopologyFeature*,const int&);
    friend void removeWalkerCluster_(TopologyFeature *,const int&);
  };


  }

#endif  // MYTHICAL_CLUSTER_HPP
