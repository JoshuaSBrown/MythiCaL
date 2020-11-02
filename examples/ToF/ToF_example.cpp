
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

int main() {
  
  const double nm_to_m = 1E-9; // m/nm

  // Define the dimensions of our system in terms of lattice sites
  const int len = 50;
  const int wid = 80;
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

  // Calculate rates between neighboring sites
  double cutoff_dist = 2.0; // nm
  double lambda = 0.01; // eV
  double Temp = 300; // K 

  // Here we are have H_AB = A * exp( -alpha * r_ij ) 
  double alpha = 6; // 1/nm
  double A = 8; // eV
  auto marcus = myct::Marcus(lambda, Temp);
 
  double electric_field = 0.1; // eV/nm

  unordered_map<int,unordered_map<int,double>> distances = lattice.getNeighborDistances(cutoff_dist);

  int debug_site = 0;
  // Here we are now going to assign rates based on the Marcus rate equation
  unordered_map<int,unordered_map<int,double>> rates;
  for ( auto site_neighs : distances ) {
    int site_i = site_neighs.first;
    int x_i = lattice.getX(site_i);
    double energy_i = site_energies.at(site_i);
    for ( auto site_dist : site_neighs.second ) {
      int site_j = site_dist.first;
      int x_j = lattice.getX(site_j);
      double x_diff_dist = static_cast<double>(x_j - x_i) * inter_site_dist; 
      double energy_j = site_energies.at(site_j) * electric_field * x_diff_dist;
      double r_ij = site_dist.second;
      double H_AB = A*std::exp( -1.0 * alpha * r_ij ); 
      double rate_ij = marcus.getRate(energy_i, energy_j, H_AB);
      double rate_ji = marcus.getRate(energy_j, energy_i, H_AB);
      rates[site_i][site_j] = rate_ij;
      rates[site_j][site_i] = rate_ji;
      if ( site_i > debug_site ) debug_site = site_i;
      if ( site_j > debug_site ) debug_site = site_j;
    }
  }

  std::cout << "Max site id in rates " << debug_site << std::endl;
  // Now we are randomly going to pick sites on the first plane of our lattice
  int num_charges = 100; 
  unordered_set<int> siteOccupied;
  while ( siteOccupied.size() < num_charges ) {
    int site = lattice.getRandomSite(myct::Cubic::Plane::YZ, 0);
    siteOccupied.insert(site);
  }
  
  class Hole : public my::Walker {};

  // Now we are going to place charges on the randomly picked sites
  vector<pair<int,shared_ptr<my::Walker>>> holes;
  holes.reserve(num_charges); 
  int charge_id = 0;
  for ( int site_id : siteOccupied ) {
    holes.emplace_back(charge_id,shared_ptr<my::Walker>(new Hole));
    holes.back().second->occupySite(site_id);
    ++charge_id;
  }

  int data_samples = 200; // The number of times we want to measure the current 
  double cutoff_time = 1E-7; // seconds 
  double current_time_sample_increment = cutoff_time/static_cast<double>(data_samples);
  double sample_time = current_time_sample_increment;
  
  int random_number_seed = 1943;
  my::CoarseGrainSystem CGsystem;
  CGsystem.setRandomSeed(random_number_seed);
  CGsystem.setTimeResolution(sample_time);
  CGsystem.initializeSystem(rates);
  CGsystem.initializeWalkers(holes);

  vector<double> transient_current(data_samples,0.0);
  // Calculate Walker dwell times and sort 
  my::Queue walker_global_times;
  {

    mt19937 random_number_generator;
    random_number_generator.seed(3);
    uniform_real_distribution<double> distribution(0.0,1.0);

    std::cout << "Dwell times" << std::endl;
    for(int walker_index=0; walker_index < num_charges; ++walker_index){
      std::cout << holes.at(walker_index).second->getDwellTime() << std::endl;
      walker_global_times.add(pair<int,double>(walker_index,holes.at(walker_index).second->getDwellTime()));
    }
    walker_global_times.sort();
  }// Calculate walker dwell times and sort
  assert(walker_global_times.at(0)->second<cutoff_time);

  std::cout << std::endl;
  std::cout << "Transient Current" << std::endl;

  int current_index = 0; 
  while(!walker_global_times.size() && walker_global_times.begin()->second<cutoff_time){
    double deltaX = 0.0;
    while(!walker_global_times.empty() && walker_global_times.begin()->second<sample_time){
      auto walker_index = walker_global_times.begin()->first;
      std::shared_ptr<my::Walker>& hole = holes.at(walker_index).second;
      int siteId = hole->getIdOfSiteCurrentlyOccupying();
      int old_x_pos = lattice.getX(siteId);
      CGsystem.hop(walker_index,hole);
      siteId = hole->getIdOfSiteCurrentlyOccupying();
      int new_x_pos = lattice.getX(siteId);
      deltaX+=static_cast<double>(new_x_pos-old_x_pos);
//      std::cout << "new x pos " << new_x_pos << " delta X " << deltaX <<  std::endl;
      // Update the dwell time
      walker_global_times.begin()->second += hole->getDwellTime();
      // reorder the walkers based on which one will move next
      if(new_x_pos == lattice.getLength()-1){
        std::cout << "Removing walker" << walker_index << std::endl;
        CGsystem.removeWalkerFromSystem(walker_index,hole);
        walker_global_times.pop_front();
      }
      walker_global_times.sort(compareSecondItemOfPair);
    }
    transient_current.at(current_index) = deltaX*nm_to_m/current_time_sample_increment;
    std::cout << sample_time << " " << transient_current.at(current_index) << std::endl;
    ++current_index;
    sample_time+=current_time_sample_increment;
  }

  // Write to file
  std::ofstream fid;
  fid.open("transient_current.txt");
  if( !fid.is_open() ){
    throw std::runtime_error("Error opening transient_current.txt file to write too");
  } 
  fid << "Transient Current" << std::endl;
  fid << transient_current.size() << std::endl;
  sample_time = 0.0;
  for( auto current : transient_current){
    sample_time+=current_time_sample_increment;
    fid << sample_time << " " << current << std::endl;
  }  
  fid.close();
  return 0;
}
