#include <cassert>
#include <cmath>
#include <iostream>

#include "../../../include/kmccoursegrain/kmc_constants.hpp"
#include "../../../include/kmccoursegrain/kmc_particle.hpp"

using namespace std;
using namespace kmccoursegrain;

int main(void) {

  cout << "Testing: Particle constructor" << endl;
  { KMC_Particle particle; }

  cout << "Testing: Particle occupySite" << endl;
  {
    KMC_Particle particle;
    particle.occupySite(0);
  }

  cout << "Testing: Particle getIdOfSiteCurrentlyOccupying" << endl;
  {
    KMC_Particle particle;

    particle.occupySite(2);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 2);
    particle.occupySite(2);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 2);
    particle.occupySite(4);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 4);
    particle.occupySite(1);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 1);
    particle.occupySite(0);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 0);
    particle.occupySite(3);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 3);
    particle.occupySite(4);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 4);
    particle.occupySite(1);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 1);
    particle.occupySite(4);
    assert(particle.getIdOfSiteCurrentlyOccupying() == 4);
  }

  cout << "Testing: Particle set and get PotentialSite" << endl;
  {
    KMC_Particle particle;
    particle.setPotentialSite(4);
    assert(particle.getPotentialSite() == 4);
  }

  cout << "Testing: Particle set and get DwellTime" << endl;
  {
    KMC_Particle particle;
    particle.setDwellTime(124.0);
    assert(static_cast<int>(round(particle.getDwellTime())) == 124);
  }
  return 0;
}
