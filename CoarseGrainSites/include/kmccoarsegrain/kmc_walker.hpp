#ifndef KMCCOARSEGRAIN_KMC_WALKER_HPP
#define KMCCOARSEGRAIN_KMC_WALKER_HPP
#include "kmc_constants.hpp"
#include <list>
#include <map>
#include <vector>

namespace kmccoarsegrain {

/**
 * \brief Walker class meant to be inherited
 *
 * This class is meant to be inherited, it keeps track of the current site that
 * is occupied as well as the dwell time and the next potential site.
 **/
class KMC_Walker {
 public:
  KMC_Walker();

  /**
   * \brief Record the id of the site the walker currently resides on
   *
   * \param[in] siteId id of the site the walker is to occupy
   **/
  void occupySite(const int & siteId) { current_site_ = siteId; }

  /**
   * \brief Grabs the id of the site the walker currently resides on
   *
   * \return site id
   **/
  int getIdOfSiteCurrentlyOccupying() const;

  /**
   * \brief Record the site the walker will move to next
   *
   * \param[in] site id the walker will try to move to if it remains
   * unoccupied.
   **/
  void setPotentialSite(const int & siteId) { potential_site_ = siteId; }

  /**
   * \brief Return the id of the site the walker will attempt to move to
   **/
  int getPotentialSite() const;

  /**
   * \brief Get the dwell time
   **/
  double getDwellTime() const { return dwell_time_; }

  /**
   * \brief Set the dwell time
   **/
  void setDwellTime(const double & dwell_time) { dwell_time_ = dwell_time; }

 private:
  /// The site the walker currently resides on
  int current_site_;

  /// The next site the walker will try to move to
  int potential_site_;

  /// The length of time the walker will remain on the current site
  double dwell_time_;
};
}
#endif  // KMCCOARSEGRAIN_KMC_WALKER_HPP
