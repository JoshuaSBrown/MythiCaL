
#include <mythical/queue.hpp>
#include <mythical/walker.hpp>
#include <mythical/charge_transport/cubic_lattice.hpp>
#include <mythical/charge_transport/marcus.hpp>
#include <mythical/coarsegrainsystem.hpp> 
#include <cassert> 
#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <random>
#include <unordered_map>
#include <unordered_set>

using namespace std;

namespace my = mythical;
namespace myct = mythical::charge_transport;

typedef std::pair<int,double> element_t;
typedef std::pair<int,shared_ptr<my::Walker>> walker_t;

std::vector<double> generateSiteEnergies(const int & total_num_sites) {
  std::cout << "- Generating site energies" << std::endl;
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
  return site_energies;
}

unordered_map<int,unordered_map<int,double>> 
calculateRates(const myct::Cubic & lattice, std::vector<double> site_energies) {
  std::cout << "- Calculating rates between sites." << std::endl;
  // Calculate rates between neighboring sites
  double cutoff_dist = 2.0; // nm
  double lambda = 0.02; // eV
  double Temp = 300; // K 
 
  // Here we are have H_AB = A * exp( -alpha * r_ij ) 
  double alpha = 6; // 1/nm
  double A = 8; // eV
  auto marcus = myct::Marcus(lambda, Temp);
 
  double electric_field = 0.1; // eV/nm
  double charge = 1.0; // q because hole

  unordered_map<int,unordered_map<int,double>> distances = lattice.getNeighborDistances(cutoff_dist);
  // Here we are now going to assign rates based on the Marcus rate equation
  unordered_map<int,unordered_map<int,double>> rates;
  for ( auto site_neighs : distances ) {
    int site_i = site_neighs.first;
    int x_i = lattice.getX(site_i);
    double energy_i = site_energies.at(site_i);
    for ( auto site_dist : site_neighs.second ) {
      int site_j = site_dist.first;
      int x_j = lattice.getX(site_j);
      double x_diff_dist = static_cast<double>(x_j - x_i) * lattice.getLatticeSpacing(); 
      double energy_j = site_energies.at(site_j) - electric_field * x_diff_dist * charge;
      double r_ij = site_dist.second;
      double H_AB = A*std::exp( -1.0 * alpha * r_ij ); 
      double rate_ij = marcus.getRate(energy_i, energy_j, H_AB);
      double rate_ji = marcus.getRate(energy_j, energy_i, H_AB);
      rates[site_i][site_j] = rate_ij;
      rates[site_j][site_i] = rate_ji;
    }
  }
  return rates;
}

class Hole : public my::Walker {};

vector<walker_t> populateLattice(int num_charges, myct::Cubic & lattice ) {
  std::cout << "- Populating first plane of lattice with charges" << std::endl;
  unordered_set<int> siteOccupied;
  while ( siteOccupied.size() < num_charges ) {
    int site = lattice.getRandomSite(myct::Cubic::Plane::YZ, 0);
    siteOccupied.insert(site);
  }

  // Now we are going to place charges on the randomly picked sites
  vector<walker_t> holes;
  holes.reserve(num_charges); 
  int charge_id = 0;
  for ( int site_id : siteOccupied ) {
    holes.emplace_back(charge_id,shared_ptr<my::Walker>(new Hole));
    holes.back().second->occupySite(site_id);
    ++charge_id;
  }
  return holes;
}

my::Queue createQueue(const vector<walker_t> & holes, double cutoff_time) {
  std::cout << "- Creating queue for holes." << std::endl;
  my::Queue walker_global_times;
  mt19937 random_number_generator;
  random_number_generator.seed(4);
  uniform_real_distribution<double> distribution(0.0,1.0);

  for(int walker_index=0; walker_index < holes.size(); ++walker_index){
    walker_global_times.add(pair<int,double>(walker_index,holes.at(walker_index).second->getDwellTime()));
  }
  walker_global_times.sort();
  assert(walker_global_times.at(0).second<cutoff_time);

  return walker_global_times;
}

void printTransientCurrent(
    const std::vector<double> & transient_current, 
    const std::vector<double> & charges_left, 
    const double current_time_sample_increment)
{
  std::cout << "- Writing transient current to file." << std::endl;
  // Write to file
  std::ofstream fid;
  fid.open("transient_current.txt");
  if( !fid.is_open() ){
    throw std::runtime_error("Error opening transient_current.txt file to write too");
  } 
  fid << "Transient Current" << std::endl;
  fid << transient_current.size() << std::endl;
  double sample_time = 0.0;
  for( int ind = 0; ind < transient_current.size(); ++ind){
    sample_time += current_time_sample_increment;
    fid << sample_time << " " << transient_current.at(ind) << " " << charges_left.at(ind) << std::endl;
  }  
  fid.close();
}

int main() {
  const double nm_to_cm = 1E-7; // cm/nm

  // Define the dimensions of our system in terms of lattice sites
  const int len = 200;
  const int wid = 80;
  const int hei = 80;
  const double inter_site_dist = 1; // nm
  const int total_num_sites = len * wid * hei; 
  myct::BoundarySetting x_b = myct::BoundarySetting::Fixed;
  myct::BoundarySetting y_b = myct::BoundarySetting::Periodic;
  myct::BoundarySetting z_b = myct::BoundarySetting::Periodic;
  myct::Cubic lattice(len, wid, hei, inter_site_dist, x_b, y_b, z_b);

  std::vector<double> site_energies = generateSiteEnergies(total_num_sites);

  auto rates = calculateRates(lattice,site_energies);

  // Now we are randomly going to pick sites on the first plane of our lattice
  int num_charges = 200; 
  auto holes = populateLattice(num_charges, lattice);

//  int data_samples = 400; // The number of times we want to measure the current 
  double cutoff_time = 0.4E-8; // seconds 
  double current_time_sample_increment = 2E-9/static_cast<double>(300);
  double sample_time = current_time_sample_increment;
  
  int random_number_seed = 1943;
  my::CoarseGrainSystem CGsystem;
  CGsystem.setRandomSeed(random_number_seed);
  CGsystem.setTimeResolution(sample_time);
  CGsystem.initializeSystem(rates);
  CGsystem.initializeWalkers(holes);

  my::Queue walker_global_times = createQueue(holes, cutoff_time);
  vector<double> transient_current;
  vector<double> charges_left;
  transient_current.reserve(300*2);
  charges_left.reserve(300*2);
  // Calculate Walker dwell times and sort 

  std::cout << "- Beginning loop." << std::endl;
  int current_index = 0; 
  int charges_remain = num_charges; 
  while(walker_global_times.size() && walker_global_times.at(0).second<cutoff_time){
    double deltaX = 0.0;
    while(walker_global_times.size() && walker_global_times.at(0).second < sample_time){
      element_t walker_time = walker_global_times.pop_current();
      std::shared_ptr<my::Walker>& hole = holes.at(walker_time.first).second;
      int siteId = hole->getIdOfSiteCurrentlyOccupying();
      int old_x_pos = lattice.getX(siteId);
      CGsystem.hop(walker_time.first,hole);
      siteId = hole->getIdOfSiteCurrentlyOccupying();
      int new_x_pos = lattice.getX(siteId);
      deltaX+=static_cast<double>(new_x_pos-old_x_pos);
      // Update the dwell time
      walker_time.second += hole->getDwellTime();
      // reorder the walkers based on which one will move next
      if(new_x_pos < lattice.getLength()-1){
        walker_global_times.sortedAdd(walker_time);
      } else {
        std::cout << "- Removing walker " << walker_time.first << " time " << walker_time.second << std::endl;
        CGsystem.removeWalkerFromSystem(walker_time.first,hole);
        --charges_remain;
      }
    }
    std::cout << "DeltaX " << deltaX << std::endl;
    std::cout << "num charges " << num_charges << std::endl;
    // Current Density velocity of charges * current_density [ J/cm^2 ]
    double q = 1.602176634E-19;
    double velocity = deltaX*nm_to_cm / sample_time; 
    std::cout << "Velocity " << velocity << std::endl;
    double num_charges = walker_global_times.size();
    double volume = static_cast<double>(len*wid*hei) * nm_to_cm*nm_to_cm*nm_to_cm;
    double current_density = q * static_cast<double>(num_charges) / volume; 
    double transient_cur = current_density * velocity;
    std::cout << "Current density " << transient_cur << std::endl;
    transient_current.push_back(transient_cur);
    charges_left.push_back(charges_remain);
    std::cout << sample_time << " " << transient_current.at(current_index) << std::endl;
    ++current_index;
    sample_time+=current_time_sample_increment;
  }
  std::cout << "- Charges remaining " << walker_global_times.size() << std::endl;
  printTransientCurrent(transient_current,charges_left, current_time_sample_increment);

  return 0;
}
