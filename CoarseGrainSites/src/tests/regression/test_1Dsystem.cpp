#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <string>
#include <random>
#include <map>
#include <unordered_map>
#include <set>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "../../../include/kmccoarsegrain/kmc_coarsegrainsystem.hpp"
#include "../../../include/kmccoarsegrain/kmc_walker.hpp"

using namespace std;
using namespace std::chrono;
using namespace kmccoarsegrain;

bool compareSecondItemOfPair(const pair<int,double> &x, const pair<int,double> & y){
  return x.second<y.second;
}

int main(int argc, char* argv[]){

  if(argc!=7){
    cerr << "To run the program correctly you must provide the " << endl;
    cerr << "following parameters: " << endl;
    cerr << endl;
    cerr << "distance   - integer defines the length of the 1d system" << endl;
    cerr << "             in terms of the number of sites. [nm] " << endl;
    cerr << "seed       - integer value defines the random number seed" << endl;
    cerr << "threshold  - integer value defines minimum threshold of " << endl;
    cerr << "             how often the simulation will try to coarse " << endl;
    cerr << "             grain." << endl;
    cerr << "time       - a double defining how long the simulation " << endl;
    cerr << "             will run for." << endl;
    cerr << "sample     - how often the current is sampled with respect"<< endl; 
    cerr << "rate         to the cutoff time." <<endl;
    cerr << "field      - the electric field across the system. [V/cm]" << endl;
    cerr << endl;
    cerr << "To run:" << endl;
    cerr << endl;
    cerr << "./performance_test_crude_vs_coarsegraining distance seed threshold time sampe_rate field" << endl;
    cerr << endl;
    return -1;
  }

  int distance = stoi(string(argv[1]));
  int seed     = stoi(string(argv[2]));
  int threshold = stoi(string(argv[3]));
  double cutoff_time = stod(string(argv[4]));
  int sample_rate = stoi(string(argv[5]));
  double field = stod(string(argv[6]));

  cout << endl;
  cout << "Parameters passed in:" << endl;
  cout << endl;
  cout << "distance:    " << distance << endl;
  cout << "seed:        " << seed << endl;
  cout << "threshold:   " << threshold << endl;
  cout << "time:        " << cutoff_time << endl;
  cout << "sample rate: " << sample_rate << endl;
  cout << "field:       " << field << endl;
  cout << endl;

  double field_nm = field*10E-7; // eV/nm
  double nm_to_m = 1E-9; // nm/m

  cout << "field [eV/nm]: " << field_nm << endl;

  double reorganization_energy = 0.01;
  // Record setup time
  int number_of_trap_sites = distance/10;
  cout << "Total number of trap sites " << number_of_trap_sites << endl;
  vector<double> energies(distance,0.0);
  {
    // Setup the energies so that two low energy sites appear every 2 sites
    for(int index=0;index<distance;++index){
      if(index%8==0 || index%9==1){
        energies.at(index)=-reorganization_energy*5;
      }
    }
  }

  unordered_map<int,unordered_map<int,double>> rates;
  unordered_map<int,vector<int>> neighbors;
  {

    double J = 0.01;
    double kBT = 0.025;
    cout << "Calculating rates using Semiclassical Marcus theory assuming: " << endl;
    cout << endl;
    cout << "reoganization energy lambda:          " << reorganization_energy << endl;
    cout << "transfer integral J:                  " << J << endl;
    cout << "Boltzmann constant * temperature kBT: " << kBT << endl;

    double hbar = pow(6.582,-16);
    double pi = 3.14;

    // Define marcus coefficient
    double coef = 2*pi/hbar*pow(J,2.0)*1/pow(4*pi*kBT,1.0/2.0);

    for(int x=0; x<distance; ++x){
      // Define neighbors
      int xlow = x;
      int xhigh = x;

      if(xlow-1>0) --xlow;
      if(xhigh+1<distance) ++xhigh;

      int siteId = x; 
      for( int x2 = xlow; x2<=xhigh; ++x2){

        double xdiff = static_cast<double>(x2-x);
        // Sign change is because a positive field should lower the energy
        double field_energy = -1.0*xdiff*field_nm;
        assert(x2>=0);
        assert(x2<distance);
        int neighId = x2;
        if(siteId!=neighId){
          neighbors[siteId].push_back(neighId);
          double deltaE = energies.at(neighId)-energies.at(siteId)-field_energy;
          double exponent = -pow(reorganization_energy-deltaE,2.0)/(4.0*reorganization_energy*kBT);
          rates[siteId][neighId] = coef*exp(exponent);
        }
      }
      assert(rates[siteId].size()!=0);          
    }
  }

  // Run coarse grained Monte Carlo
  cout << "Running coarse grained Monte Carlo" << endl;
  {
    // greating map with pointer to rates
    unordered_map< int, unordered_map< int, double *>> rates_to_neighbors;
    {
      for(auto site_rates : rates){
        for( auto neigh_rate : site_rates.second){
          rates_to_neighbors[site_rates.first][neigh_rate.first] =&(rates[site_rates.first][neigh_rate.first]);
        }
      }
    }
    // Run the coarse grain simulation
    {
      double current_time_sample_increment = cutoff_time/static_cast<double>(sample_rate);
      double sample_time = current_time_sample_increment;

      KMC_CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(seed);
      CGsystem.setMinCoarseGrainIterationThreshold(threshold);
      CGsystem.setTimeResolution(sample_time);
      CGsystem.initializeSystem(rates_to_neighbors);

      class Electron : public KMC_Walker {};
      // Only a single electron is created but the simulation is repeated 
      // several times
      int repetitions = 8;
      vector<vector<double>> all_transient_currents;
      cout << "Repetitions " << repetitions << endl;

      for(int rep=0; rep<repetitions;++rep){
        // Create the electrons using the KMC_Walker class
        vector<KMC_Walker> electrons;        
        Electron electron;
        electrons.push_back(electron);
        electrons.at(0).occupySite(0);
        CGsystem.initializeWalkers(electrons);

        vector<double> transient_current(sample_rate,0.0);
        // Calculate Walker dwell times and sort 
        list<pair<int,double>> walker_global_times;
        {

          mt19937 random_number_generator;
          random_number_generator.seed(3);
          uniform_real_distribution<double> distribution(0.0,1.0);

          walker_global_times.push_back(pair<int,double>(0,electrons.at(0).getDwellTime()));
          walker_global_times.sort(compareSecondItemOfPair);
        }// Calculate walker dwell times and sort
        assert(walker_global_times.begin()->second<cutoff_time);

        int current_index = 0; 
        while(walker_global_times.size()>0 && walker_global_times.begin()->second<cutoff_time){
          double deltaX = 0.0;
          while(walker_global_times.size()>0 && walker_global_times.begin()->second<sample_time){

            auto walker_index = walker_global_times.begin()->first;
            KMC_Walker& electron = electrons.at(walker_index); 
            int siteId = electron.getIdOfSiteCurrentlyOccupying();
            int old_x_pos = siteId;
            CGsystem.hop(electron);
            siteId = electron.getIdOfSiteCurrentlyOccupying();
            int new_x_pos = siteId;
            deltaX+=static_cast<double>(new_x_pos-old_x_pos);
            // Update the dwell time
            walker_global_times.begin()->second += electron.getDwellTime();
            // reorder the walkers based on which one will move next
            if(new_x_pos==(distance-1)){
              cout << "a. Repetition " << rep << " final position " << electrons.at(0).getIdOfSiteCurrentlyOccupying() << endl;
              CGsystem.removeWalkerFromSystem(electron);
              walker_global_times.pop_front();
            }
          }

          transient_current.at(current_index) = deltaX*nm_to_m/current_time_sample_increment;
          ++current_index;
          sample_time+=current_time_sample_increment;

        }

        if(walker_global_times.size()!=0){
          cout << "b. Repetition " << rep << " final position " << electrons.at(0).getIdOfSiteCurrentlyOccupying() << endl;
          CGsystem.removeWalkerFromSystem(electrons.at(0));
          walker_global_times.pop_front();
        }
        
        all_transient_currents.push_back(transient_current);
        cout << "end of for loop" << endl;
      } // End of the reps 

      auto clusters = CGsystem.getClusters();
      cout << "Number of clusters found " << clusters.size() << endl;

      // Calulate the mean of the transient currents
      double sum = 0.0;
      double count = 0.0;
      for( auto trans : all_transient_currents){
        for(auto current : trans){
          count+=1.0;
          sum+=current;
        }
      }
      double mean = sum/count;

      // Caluculate the standard deviation
      double sum_squares=0.0;
      for( auto trans : all_transient_currents){
        for(auto current : trans){
          sum_squares+=pow(current-mean,2.0);
        }
      }
      double standard_deviation = pow(sum_squares/count,0.5);
      cout << "Standard deviation " << standard_deviation << endl;
      cout << "Mean " << mean << endl;

    }// End of the Coarse grain simulation 
    

  } // End of coarse grain Monte Carlo

  return 0;
}
