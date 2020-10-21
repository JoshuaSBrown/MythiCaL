#include <cassert>
#include <iostream>

#include "mythical/charge_transport/cubic_lattice.hpp"

using namespace std;
using namespace mythical;
using namespace mythical::lattice;

int main(void) {

  cout << "Testing: Cubic Lattice Constructor 1" << endl;
  { Cubic lattice; }

  cout << "Testing: Cubic Lattice Constructor 2" << endl;
  { Cubic lattice(1, 1, 2); }

  cout << "Testing: Cubic Lattice Constructor 3" << endl;
  { Cubic lattice(2, 3, 4, 2.0); }

  return 0;
}
