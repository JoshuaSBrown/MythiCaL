#ifndef KMCCOARSEGRAIN_KMC_SITE_HPP
#define KMCCOARSEGRAIN_KMC_SITE_HPP
#include <iostream>
#include <unordered_map>
#include <math.h>
#include <memory>
#include <random>
#include <vector>

#include "kmc_topology_feature.hpp"

namespace kmccoarsegrain {

/**
 * \brief Site Class
 *
 * This class keeps track of all information related to a site and it's
 * neighbors. It is an internal class meaning it is not meant to be used by
 * the public. It does not store rates to the neighboring sites locally, this
 * is to make it possible for an end user to alter the rates externally
 * without having to touch this class.
 **/
class KMC_Site : public KMC_TopologyFeature {
 public:
  KMC_Site();

  ~KMC_Site();
  /**
   * \brief Sets the rates to sites neighboring this site
   *
   * This function does not remove any of the previously stored rates
   * however, it may very well overwrite a rate to a previously stored
   * neighboring site.
   *
   * \param[in] neighRates Stores the site id of the neighbor with a pointer
   * to the rate going to the neighboring site.
   **/
  void setRatesToNeighbors(std::unordered_map<int, double>& neighRates);

  /**
   * \brief Add a rate to a neighboring site
   *
   * Will add a rate to a neighboring site. Note that an error will be
   * thrown if a rate has already been added for that particular neighbor.
   * In such a case you should call resetNeighRate instead.
   *
   * \param[in] neighborRate The first int is the site id of neighboring
   * site, this is followed by a pointer to actual rate.
   **/
  void addNeighRate(const std::pair<int, double*> neighRate);

  /**
   * \brief Reset the rate to a neighboring site
   *
   * In the case that you want to change a rate or add a new one you can
   * call reset rate.
   *
   * \param[in] neighRate the first int is the id of the neighboring site,
   * the double pointer will point to the rate.
   **/
  void resetNeighRate(const std::pair<int, double*> neighRate);

  /**
   * \brief Is the site a neighbor
   *
   * Determines if the site is currently considered a neighbor.
   *
   * \param[in] neighSiteId is the id of the thought to be neighbor
   *
   * \return True if it is a neighbor and False if it is not
   **/
  bool isNeighbor(const int neighSiteId) const {
    return neighRates_.count(neighSiteId);
  }

  /**
   * \brief Determine if the site is part of a cluster
   *
   * \return True if is currently part of a cluster, False otherwise
   **/
  bool partOfCluster() const { return cluster_id_ > -1; }

  /**
   * \brief Get all the rates to the neighboring sites
   *
   * \return Returns the rates to all the neighboring sites as doubles
   **/
  std::vector<double> getRateToNeighbors() const;

  /**
   * \brief Find the rate to a neighboring site
   *
   * Return the rate to a neighboring site if it exists. If the supposed
   * site is not actually a neighbor will return an error.
   *
   * \param[in] neighSiteId the id of the site thought to be a neighbor
   *
   * \return The rate from the site to the neighbor
   **/
  double getRateToNeighbor(const int & neighSiteId) const;
  double getFastestRate();

  /**
   * \brief Get the ids of all the neighboring sites
   *
   * \return A vector of all the neighboring site ids
   **/
  std::vector<int> getNeighborSiteIds() const;

  /**
   * \brief Returns the id of a neighboring site
   *
   * Again this is based on a random number generator. A neighboring site id
   * is returned based on the rates to the neighboring site. If the rate to
   * the neighbor is large than it will be more likely that the site id
   * associated with that neighbor will be returned. Will throw an error if
   * the site has no neighbors.
   *
   * \return site id of a neigboring site
   **/
  int pickNewSiteId(const int & ) override;
  int pickNewSiteId() override;

  /**
   * \brief Return the id of the cluster the site is attached too
   *
   * \return the id of the cluster, a value of -1 means it is not attached
   * to any cluster
   **/
  int getClusterId() const { return cluster_id_; }

  /**
   * \brief Sets the id of the cluster that the site is attached to
   *
   * \param[in] clusterId the id of the cluster
   **/
  void setClusterId(const int cluster_id) { cluster_id_ = cluster_id; }

  /**
   * \brief Returns the probability of hopping to a neighboring site
   *
   * The higher the rate to the neighboring site the larger the probability
   * should be.
   *
   * \param[in] neighSiteId this is the id of the neighboring site
   *
   * \return the probabiliy of hopping to the neighbor
   **/
  double getProbabilityOfHoppingToNeighboringSite(const int & neighSiteId);

  std::unordered_map<int,const double *> & getNeighborsAndRates();
  const std::unordered_map<int,const double *>& getNeighborsAndRatesConst() const;

  /**
   * \brief Gets the ids and the probabilities of hopping to neighbors
   *
   * \return A map the first int is the neighbor site id, the double is the
   * probability that a charge located on the site will hop to one of these
   * neighbors
   **/
  std::vector<std::pair<int, double>> getProbabilitiesAndIdsOfNeighbors() const;

  /**
   * \brief Prints the output of the site
   **/
  friend std::ostream& operator<<(std::ostream& os,
                                  const kmccoarsegrain::KMC_Site& site);

   private:
  /**
   * \brief Contains the probability of hopping to each neighbor
   **/
  std::vector<std::pair<int, double>> probabilityHopToNeighbor_;

  /**
   * \brief Stores pointers to the rates to each of the neighboring sites
   **/
  std::unordered_map<int,const double*> neighRates_;

  /**
   * \brief Stores the id of the cluster the site is a part of
   **/
  int cluster_id_;

  /**
   * \brief Distribution to be used when picking random numbers
   **/
  std::uniform_real_distribution<double> randomDistribution_;

  /**
   * \brief This function calculates the probability of hopping to each
   * neighboring site
   *
   * The probabilities are stored in the probHopToNeighbor_ map
   **/
  void calculateProbabilityHopToNeighbors_();

  /**
   * \brief Calculates the escapeTimeConstant_
   **/
  void calculateDwellTimeConstant_();

  double getSumOfRates_();

};

}

#endif  // KMCCOARSEGRAIN_KMC_SITE_HPP
