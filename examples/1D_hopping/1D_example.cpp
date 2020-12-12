#include <mythical/coarsegrainsystem.hpp>
#include <mythical/walker.hpp>
    
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

namespace my = mythical;

typedef std::pair<int,shared_ptr<my::Walker>> walker_t;

unordered_map<int,unordered_map<int, double>> 
generateRates(const int number_of_sites) {
  unordered_map<int,unordered_map<int, double>> rates;
  for(int i=1 ; i<number_of_sites; ++i){
    rates[i][i+1]=3.0;
    rates[i+1][i]=1.0;
  }
  return rates;
}

int main(){

  int number_of_sites = 20;

  unordered_map<int,unordered_map<int, double>> rates =  generateRates(number_of_sites);

  my::CoarseGrainSystem CGsystem;
  CGsystem.setTimeResolution(3.0);
  CGsystem.initializeSystem(rates);

  class Charge : public my::Walker {};

  int charge_id;
  vector<walker_t> charges;
  charges.emplace_back(charge_id, std::shared_ptr<my::Walker>( new Charge));
  charges.back().second->occupySite(1);
  CGsystem.initializeWalkers(charges);

  const int iterations = 50;

  ofstream traj_file;
  traj_file.open("Trajectory.txt");
  traj_file << "Number of Sites: " << number_of_sites << endl;
  traj_file << "Number of iterations: " << iterations << endl;
  traj_file << endl;
  traj_file << "Site    Dwell Time" << std::endl;

  for(int iteration = 0; iteration<iterations;++iteration){
    CGsystem.hop(charges.at(0));
    const int site_id = charges.at(0).second->getIdOfSiteCurrentlyOccupying();
    const double dwell_time = charges.at(0).second->getDwellTime();
    cout << "Charge occupying site: " << site_id << " dwell time " << dwell_time << endl;
    traj_file << site_id << "      " << dwell_time << endl;
  }

  CGsystem.removeWalkerFromSystem(charges.at(0));
  traj_file.close();

  return 0;
}
