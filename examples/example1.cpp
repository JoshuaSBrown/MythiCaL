
#include <mythical/charge_transport/cubic_lattice.hpp>

#include <cassert> 
#include <vector>
#include <random>
#include <unordered_map>
#include <unordered_set>

using namespace std;

namespace myct = mythical::charge_transport;

int main() {

  // Define the dimensions of our system in terms of lattice sites
  const int len = 100;
  const int wid = 50;
  const int hei = 50;
  const double inter_site_dist = 1.0;
  const int total_num_sites = len * wid * hei; 
  myct::BoundarySetting x_b = myct::BoundarySetting::Fixed;
  myct::BoundarySetting y_b = myct::BoundarySetting::Periodic;
  myct::BoundarySetting z_b = myct::BoundarySetting::Periodic;
  myct::Cubic lattice(len, wid, hei, inter_site_dist, x_b, y_b, z_b);

  // Define our DOS as a normal distribution 
  const double mean = -1.2;
  const double std_deviation = 0.07;
  std::normal_distribution<double> distribution(mean, std_deviation);
 
  // Randomly assign energies from our DOS to our lattice sites
  std::default_random_engine generator;
  std::vector<double> site_energies(total_num_sites);
  for(size_t index = 0; index < total_num_sites; ++index){
    site_energies.at(index) = distribution(generator);
  }

  unordered_map<int,unordered_map<int,double>> rates;
  unordered_map<int,vector<int>> neighbors;

  unordered_set<int> siteOccupied;

  class Hole : public Walker {};
  
  vector<pair<int,Holes>> holes;        

  return 0;
}
