#include "../../include/kmccoursegrain/kmc_particle.hpp"
#include <cassert>
#include <stdexcept>

using namespace std;

namespace kmccoursegrain {

/****************************************************************************
 * External Methods
 ****************************************************************************/
KMC_Particle::KMC_Particle() {
  current_site_ = constants::unassignedId;
  potential_site_ = constants::unassignedId;
  dwell_time_ = -1.0;
}

int KMC_Particle::getIdOfSiteCurrentlyOccupying() const {
  if (current_site_ == constants::unassignedId) {
    throw runtime_error(
        "Cannot get current site as it has not yet been "
        "assigned. You many need to first initialize the particle");
  }
  return current_site_;
}

int KMC_Particle::getPotentialSite() const {
  if (potential_site_ == constants::unassignedId) {
    throw runtime_error(
        "Cannot get potential site as it has not yet been "
        "assigned. You many need to first initialize the particle");
  }
  return potential_site_;
}
}
