#ifndef KMCCOURSEGRAIN_KMC_PARTICLE_HPP
#define KMCCOURSEGRAIN_KMC_PARTICLE_HPP

#include "kmc_constants.hpp"
#include <list>
#include <map>
#include <vector>

namespace kmccoursegrain {

/**
 * \brief Particle class meant to be inherited
 *
 * This class is meant to be inherited, it keeps track of the current site that
 * is occupied as well as the dwell time and the next potential site.
 **/
class KMC_Particle {
 public:
  KMC_Particle();

  /**
   * \brief Record the id of the site the particle currently resides on
   *
   * \param[in] siteId id of the site the particle is to occupy
   **/
  void occupySite(int siteId) { current_site_ = siteId; }

  /**
   * \brief Grabs the id of the site the particle currently resides on
   *
   * \return site id
   **/
  int getIdOfSiteCurrentlyOccupying() const;

  /**
   * \brief Record the site the particle will move to next
   *
   * \param[in] site id the particle will try to move to if it remains
   * unoccupied.
   **/
  void setPotentialSite(const int siteId) { potential_site_ = siteId; }

  /**
   * \brief Return the id of the site the particle will attempt to move to
   **/
  int getPotentialSite() const;

  /**
   * \brief Get the dwell time
   **/
  double getDwellTime() const { return dwelltime_; }

  /**
   * \brief Set the dwell time
   **/
  void setDwellTime(const double dwelltime) { dwelltime_ = dwelltime; }

 private:
  /// The site the particle currently resides on
  int current_site_;

  /// The next site the particle will try to move to
  int potential_site_;

  /// The length of time the particle will remain on the current site
  double dwelltime_;
};
}
#endif  // KMCCOURSEGRAIN_KMC_PARTICLE_HPP
