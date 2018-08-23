#include <iostream>
#include <cassert>

#include <kmccoursegrain/particle.hpp>

using namespace std;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: Particle constructor" << endl;
  {
    Particle particle;
  }

  cout << "Testing: Particle getMemoryCapacity" << endl;
  {
    Particle particle;
    assert(particle.getMemoryCapacity()==0);
  }
	
  cout << "Testing: Particle setMemoryCapacity" << endl;
  {
    Particle particle;
    particle.setMemoryCapacity(5);
    assert(particle.getMemoryCapacity()==5);
    
    bool fail = false;
    try {
      particle.setMemoryCapacity(2);
    }catch(...){
      fail = true;
    }
    assert(fail);
  }

  cout << "Testing: Partilce occupySite" << endl;
  {
    Particle particle;
    particle.setMemoryCapacity(5);
    particle.occupySite(0);
  }

  cout << "Testing: Particle getVisitationFrequencyOfCurrentlyOccupiedSite" << endl;
  {
    Particle particle;
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

	return 0;
}
