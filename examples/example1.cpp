
#include <mythical/charge_transport/cubic_lattice.hpp>

#include <cassert> 
#include <vector>
#include <random>

using namespace std;

namespace myla = mythical::lattice;

int main() {

  const int len = 100;
  const int wid = 50;
  const int hei = 50;
  const double inter_site_dist = 1.0;
  myla::BoundarySetting x_b = myla::BoundarySetting::Fixed;
  myla::BoundarySetting y_b = myla::BoundarySetting::Periodic;
  myla::BoundarySetting z_b = myla::BoundarySetting::Periodic;

  myla::Cubic lattice(len, wid, hei, inter_site_dist, x_b, y_b, z_b);

  const int total_num_sites = len * wid * hei; 
  std::vector<double> site_energies(total_num_sites);

  const double mean = -1.2;
  const double std_deviation = 0.07;

  std::normal_distribution<double> distribution(mean, std_deviation);
  
  std::default_random_engine generator;
  for(size_t index = 0; index < total_num_sites; ++index){
    site_energies.at(index) = distribution(generator);
  }

  return 0;
}
