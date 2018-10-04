#include <iostream>
#include <cassert>

#include <kmccoursegrain/kmc_particle.hpp>

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
    cout << particle.getVisitationFrequencyOfCurrentlyOccupiedSite() << endl;
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

	return 0;
}
