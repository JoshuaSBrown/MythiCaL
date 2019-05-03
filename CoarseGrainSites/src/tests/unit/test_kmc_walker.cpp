#include <cassert>
#include <cmath>
#include <iostream>

#include "../../../include/kmccoarsegrain/kmc_constants.hpp"
#include "../../../include/kmccoarsegrain/kmc_walker.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void) {

  cout << "Testing: Walker constructor" << endl;
  { KMC_Walker walker; }

  cout << "Testing: Walker occupySite" << endl;
  {
    KMC_Walker walker;
    int siteId = 0;
    walker.occupySite(siteId);
  }

  cout << "Testing: Walker getIdOfSiteCurrentlyOccupying" << endl;
  {
    KMC_Walker walker;
    int siteId = 2;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 2);
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 2);
    siteId = 4;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 4);
    siteId = 1;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 1);
    siteId = 0;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 0);
    siteId = 3;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 3);
    siteId = 4;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 4);
    siteId = 1;
    walker.occupySite(siteId);
    assert(walker.getIdOfSiteCurrentlyOccupying() == 1);
    siteId = 4;
    walker.occupySite(siteId);
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
