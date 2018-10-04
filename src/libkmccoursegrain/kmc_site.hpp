#ifndef KMCCOURSEGRAIN_KMC_SITE_HPP
#define KMCCOURSEGRAIN_KMC_SITE_HPP

#include <iostream>
#include <vector>
#include <map>
#include <memory>
#include <math.h> 
#include <random>

#include "identity.hpp"

namespace kmccoursegrain {

  /**
   * \brief Site Class
   *
   * This class keeps track of all information related to a site and it's 
   * neighbors. It is an internal class meaning it is not meant to be used by
   * the public. It does not store rates to the neighboring sites locally, this
   * is to make it possible for an end user to alter the rates externally 
   * without having to touch this class. 
   **/
  class KMC_Site : public virtual Identity {

    public:	
 
      KMC_Site();

      /**
       * \brief Set the seed for the random number generator
       *
       * By default this is initialized in the constructor from the time.
       * However, having the ability to set it allows to reproducably test the
       * class.
       *
       * \param[in] seed a random number seed
       **/
      void setRandomSeed(const unsigned long seed);

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
      void setRatesToNeighbors(std::map<int const, double* > neighRates);

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
      void addNeighRate(const std::pair<const int,double * > neighRate);

      /**
       * \brief Reset the rate to a neighboring site
       *
       * In the case that you want to change a rate or add a new one you can
       * call reset rate. 
       *
       * \param[in] neighRate the first int is the id of the neighboring site,
       * the double pointer will point to the rate.
       **/
      void resetNeighRate(const std::pair<const int,double * > neighRate);

      /**
       * \brief Is the site a neighbor
       *
       * Determines if the site is currently considered a neighbor.
       *
       * \param[in] neighSiteId is the id of the thought to be neighbor
       *
       * \return True if it is a neighbor and False if it is not
       **/
      bool isNeighbor(const int neighSiteId) const 
      { return neighRates_.count(neighSiteId);}

      /**
       * \brief Determine if site is occupied by a particle
       *
       * \return True if it is occupied, False if it is not occupied
       **/
      bool siteIsOccupied() const { return siteOccupied_;}

      /**
       * \brief Determine if the site is part of a cluster
       *
       * \return True if is currently part of a cluster, False otherwise
       **/
      bool partOfCluster() const { return clusterId_>-1; }

      /**
       * \brief Get all the rates to the neighboring sites
       *
       * \return Returns the rates to all the neighboring sites as doubles
       **/
      std::vector<double> getRateToNeighbors();

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
      double getRateToNeighbor(const int neighSiteId);

      /** 
       * \brief Indicate that the site is no longer occupied by a particle
       **/
      void vacateSite(){ siteOccupied_ = false; }

      /**
       * \brief Indicate that the site is now occupied by a particle
       **/

      void occupySite(){ siteOccupied_ = true; ++totalVisitFreq_; }
      /**
       * \brief Get the ids of all the neighboring sites
       *
       * \return A vector of all the neighboring site ids
       **/
      std::vector<int> getNeighborSiteIds() const;

      /**
       * \brief Return the dwell time of the site
       *
       * \return A double which is the dwell time
       **/
      double getTimeConstant() const { return escapeTimeConstant_;}

      /**
       * \brief Return the hop time of the site
       *
       * This is different from the dwell time as it is multiplied by a random
       * number, each time this function is called it should give a different
       * time. This is where the Monte Carlo part of the simulation comes into
       * play. 
       *
       * \return A time indicating how long a particle was on the site before it
       * hopped
       **/
      double getDwellTime();

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
      int pickNewSiteId();
     
      /**
       * \brief Return the id of the cluster the site is attached too
       *
       * \return the id of the cluster, a value of -1 means it is not attached
       * to any cluster
       **/
      int getClusterId() const { return clusterId_; }

      /**
       * \brief Sets the id of the cluster that the site is attached to
       *
       * \param[in] clusterId the id of the cluster
       **/
      void setClusterId(const int clusterId) { clusterId_ = clusterId; }

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
      double getProbabilityOfHoppingToNeighboringSite(const int neighSiteId);

      /**
       * \brief Gets the ids and the probabilities of hopping to neighbors
       *
       * \return A map the first int is the neighbor site id, the double is the
       * probability that a charge located on the site will hop to one of these
       * neighbors
       **/  
      std::map<const int, double> getProbabilitiesAndIdsOfNeighbors() const;

      /**
       * \brief Prints the output of the site
       **/
      friend std::ostream& 
        operator<<(std::ostream& os, const kmccoursegrain::KMC_Site& site);
    private:

      /**
       * \brief Contains the probability of hopping to each neighbor
       **/
      std::map<int const, double> probabilityHopToNeighbor_;

      /**
       * \brief Stores pointers to the rates to each of the neighboring sites
       **/
      std::map<int const, double *> neighRates_;

      /**
       * \brief Keeps track of the total number of time the site has been
       * visited by a particle
       **/
      int totalVisitFreq_;

      /**
       * \brief Stores the id of the cluster the site is a part of
       **/
      int clusterId_;

      /**
       * \brief Keeps track of whether a site is currently occupied 
       **/
      bool siteOccupied_;

      /**
       * \brief The time constant of the site, dwell time constant 
       *
       * This value is calcualted by considering the rates off the site. If all
       * the rates off the site are very small the time constant will be large.
       * If the rates off are fast than the time constant will be small. It has
       * an inverse relation with the rates off. 
       */ 
      double escapeTimeConstant_;

      /**
       * \brief Random number seed
       **/
      unsigned long seed_;

      /**
       * \brief This is the random number engine
       *
       * This engine creates pseudo random numbers, it has been initialized 
       * from the time.
       **/
      std::mt19937 randomEngine_;

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

#endif // KMCCOURSEGRAIN_KMC_SITE_HPP
