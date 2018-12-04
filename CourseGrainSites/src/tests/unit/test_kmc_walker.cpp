#include <cassert>
#include <cmath>
#include <iostream>

#include "../../../include/kmccoursegrain/kmc_constants.hpp"
#include "../../../include/kmccoursegrain/kmc_walker.hpp"

using namespace std;
using namespace kmccoursegrain;

int main(void) {

  cout << "Testing: Walker constructor" << endl;
  { KMC_Walker walker; }

  cout << "Testing: Walker occupySite" << endl;
  {
    KMC_Walker walker;
    walker.occupySite(0);
  }

  cout << "Testing: Walker getIdOfSiteCurrentlyOccupying" << endl;
  {
    KMC_Walker walker;

    walker.occupySite(2);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 2);
    walker.occupySite(2);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 2);
    walker.occupySite(4);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 4);
    walker.occupySite(1);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 1);
    walker.occupySite(0);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 0);
    walker.occupySite(3);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 3);
    walker.occupySite(4);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 4);
    walker.occupySite(1);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 1);
    walker.occupySite(4);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 4);
  }

  cout << "Testing: Walker set and get PotentialSite" << endl;
  {
    KMC_Walker walker;
    walker.setPotentialSite(4);
    assert(walker.getPotentialSite() == 4);
  }

  cout << "Testing: Walker set and get DwellTime" << endl;
  {
    KMC_Walker walker;
    walker.setDwellTime(124.0);
    assert(static_cast<int>(round(walker.getDwellTime())) == 124);
  }
  return 0;
}
