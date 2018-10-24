#include <iostream>
#include <cmath>
#include <cassert>

#include "../../include/kmccoursegrain/kmc_particle.hpp"
#include "../../include/kmccoursegrain/kmc_constants.hpp"

using namespace std;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: Particle constructor" << endl;
  {
    KMC_Particle particle;
  }

  cout << "Testing: Particle getMemoryCapacity" << endl;
  {
    KMC_Particle particle;
    assert(particle.getMemoryCapacity()==2);
  }
	
  cout << "Testing: Particle setMemoryCapacity" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);
    assert(particle.getMemoryCapacity()==5);
    
    bool fail = false;
    try {
      particle.setMemoryCapacity(1);
    }catch(...){
      fail = true;
    }
    assert(fail);
  }

  cout << "Testing: Particle occupySite" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);
    particle.occupySite(0);
  }

  cout << "Testing: Particle getVisitationFrequencyOfCurrentlyOccupiedSite" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);     
    particle.occupySite(0);
    assert(particle.getVisitationFrequencyOfCurrentlyOccupiedSite()==1);
    particle.occupySite(1);
    assert(particle.getVisitationFrequencyOfCurrentlyOccupiedSite()==1);
    particle.occupySite(1);
    assert(particle.getVisitationFrequencyOfCurrentlyOccupiedSite()==2);

    particle.occupySite(2);
    particle.occupySite(3);
    particle.occupySite(4);
    particle.occupySite(4);
    particle.occupySite(4);
    assert(particle.getVisitationFrequencyOfCurrentlyOccupiedSite()==3);
    
  }

  cout << "Testing: Particle getMemory" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);     

    particle.occupySite(2);
    particle.occupySite(2);
    particle.occupySite(4);
    particle.occupySite(1);
    particle.occupySite(0);
    particle.occupySite(3);
    particle.occupySite(4);
    particle.occupySite(1);
    particle.occupySite(4);
    auto visits = particle.getMemory();
    
    // Site id of most recently visited to latest
    assert(visits.at(0).at(0)==4);
    assert(visits.at(1).at(0)==1);
    assert(visits.at(2).at(0)==3);
    assert(visits.at(3).at(0)==0);
    assert(visits.at(4).at(0)==2);

    // Cluster id of most recently visited site to latest
    assert(visits.at(0).at(1)==constants::unassignedId);
    assert(visits.at(1).at(1)==constants::unassignedId);
    assert(visits.at(2).at(1)==constants::unassignedId);
    assert(visits.at(3).at(1)==constants::unassignedId);
    assert(visits.at(4).at(1)==constants::unassignedId);

    // Frequency of visitation of the most recently visited site to the oldest
    assert(visits.at(0).at(2)==3); // site 4 visited 3 times
    assert(visits.at(1).at(2)==2); // site 1 visited 2 times
    assert(visits.at(2).at(2)==1); // Site 3 visited 1 time
    assert(visits.at(3).at(2)==1); // Site 0 visited 1 time
    assert(visits.at(4).at(2)==2); // Site 2 visited 2 times
    
  }

  cout << "Testing: Particle getIdOfSiteCurrentlyOccupying" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);     

    particle.occupySite(2);
    assert(particle.getIdOfSiteCurrentlyOccupying()==2);
    particle.occupySite(2);
    assert(particle.getIdOfSiteCurrentlyOccupying()==2);
    particle.occupySite(4);
    assert(particle.getIdOfSiteCurrentlyOccupying()==4);
    particle.occupySite(1);
    assert(particle.getIdOfSiteCurrentlyOccupying()==1);
    particle.occupySite(0);
    assert(particle.getIdOfSiteCurrentlyOccupying()==0);
    particle.occupySite(3);
    assert(particle.getIdOfSiteCurrentlyOccupying()==3);
    particle.occupySite(4);
    assert(particle.getIdOfSiteCurrentlyOccupying()==4);
    particle.occupySite(1);
    assert(particle.getIdOfSiteCurrentlyOccupying()==1);
    particle.occupySite(4);
    assert(particle.getIdOfSiteCurrentlyOccupying()==4);
  }

  cout << "Testing: Particle removeMemory" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);     

    particle.occupySite(2);
    particle.occupySite(2);
    particle.occupySite(4);
    particle.occupySite(1);
    particle.occupySite(0);
    particle.occupySite(3);
    particle.occupySite(4);
    particle.occupySite(1);
    particle.occupySite(4);

    particle.removeMemory(4);

    auto visits = particle.getMemory();
    
    // Site id of most recently visited to latest
    assert(visits.at(0).at(0)==1);
    assert(visits.at(1).at(0)==3);
    assert(visits.at(2).at(0)==0);
    assert(visits.at(3).at(0)==2);

    // Cluster id of most recently visited site to latest
    assert(visits.at(0).at(1)==constants::unassignedId);
    assert(visits.at(1).at(1)==constants::unassignedId);
    assert(visits.at(2).at(1)==constants::unassignedId);
    assert(visits.at(3).at(1)==constants::unassignedId);

    // Frequency of visitation of the most recently visited site to the oldest
    assert(visits.at(0).at(2)==2); // site 1 visited 2 times
    assert(visits.at(1).at(2)==1); // Site 3 visited 1 time
    assert(visits.at(2).at(2)==1); // Site 0 visited 1 time
    assert(visits.at(3).at(2)==2); // Site 2 visited 2 times
    
  }

  cout << "Testing: Particle set and get PotentialSite" << endl;
  {
    KMC_Particle particle;
    particle.setPotentialSite(4);
    assert(particle.getPotentialSite()==4);
  }

  cout << "Testing: Particle resetVisitationFrequency" << endl;
  {
    KMC_Particle particle;
    particle.setMemoryCapacity(5);     

    particle.occupySite(2);
    particle.occupySite(2);
    particle.occupySite(4);
    particle.occupySite(1);
    particle.occupySite(0);
    particle.occupySite(3);
    particle.occupySite(4);
    particle.occupySite(1);
    particle.occupySite(4);

    assert(particle.getVisitationFrequencyOfCurrentlyOccupiedSite()==3);
    assert(particle.getIdOfSiteCurrentlyOccupying()==4);
 
    particle.resetVisitationFrequency(4);
   
    assert(particle.getVisitationFrequencyOfCurrentlyOccupiedSite()==0);
    assert(particle.getIdOfSiteCurrentlyOccupying()==4);
  }

  cout << "Testing: Particle set and get DwellTime" << endl;
  {
    KMC_Particle particle;
    particle.setDwellTime(124.0);
    assert(static_cast<int>(round(particle.getDwellTime()))==124);
  }
	return 0;
}
