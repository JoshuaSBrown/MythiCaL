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

#include "../../../include/kmccoursegrain/kmc_coursegrainsystem.hpp"
#include "../../../include/kmccoursegrain/kmc_particle.hpp"

using namespace std;
using namespace std::chrono;
using namespace kmccoursegrain;

bool compareSecondItemOfPair(const pair<int,double> &x, const pair<int,double> & y){
  return x.second<y.second;
}

int main(int argc, char* argv[]){

  if(argc!=8){
    cerr << "To run the program correctly you must provide the " << endl;
    cerr << "following parameters: " << endl;
    cerr << endl;
    cerr << "distance   - integer defines the length of the 1d system" << endl;
    cerr << "             in terms of the number of sites. [nm] " << endl;
    cerr << "seed       - integer value defines the random number seed" << endl;
    cerr << "resolution - integer value defines how course the " << endl;
    cerr << "             approxiation will be." << endl;
    cerr << "threshold  - integer value defines minimum threshold of " << endl;
    cerr << "             how often the simulation will try to course " << endl;
    cerr << "             grain." << endl;
    cerr << "time       - a double defining how long the simulation " << endl;
    cerr << "             will run for." << endl;
    cerr << "sample     - how often the current is sampled with respect"<< endl; 
    cerr << "rate         to the cutoff time." <<endl;
    cerr << "field      - the electric field across the system. [V/cm]" << endl;
    cerr << endl;
    cerr << "To run:" << endl;
    cerr << endl;
    cerr << "./performance_test_crude_vs_coursegraining distance seed resolution threshold time sampe_rate field" << endl;
    cerr << endl;
    return -1;
  }

  int distance = stoi(string(argv[2]));
  int seed     = stoi(string(argv[3]));
  int resolution = stoi(string(argv[4])); 
  int particles = 1;
  int threshold = stoi(string(argv[6]));

  double cutoff_time = stod(string(argv[7]));
  int sample_rate = stoi(string(argv[8]));
  double field = stod(string(argv[9]));

  cout << endl;
  cout << "Parameters passed in:" << endl;
  cout << endl;
  cout << "distance:    " << distance << endl;
  cout << "seed:        " << seed << endl;
  cout << "resolution:  " << resolution << endl; 
  cout << "threshold:   " << threshold << endl;
  cout << "time:        " << cutoff_time << endl;
  cout << "sample rate: " << sample_rate << endl;
  cout << "field:       " << field << endl;
  cout << endl;

  double field_nm = field*10E-7; // eV/nm

  double reorganization_energy = 0.01;
  // Record setup time
  vector<double> energies(distance,0.0);
  {
    // Setup the energies so that two low energy sites appear every 2 sites
    for(int index=0;index<distance;++index){
      if(index%4==0 || index%4==1){
        energies.at(index)=-reorganization_energy*2;
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
                double field_energy = xdiff*field_nm;
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
            }
          }
          assert(rates[siteId].size()!=0);          
        }
      }
    }
  }

  // Place particles randomly on the first plane of the system for ToF simulation
  set<int> siteOccupied;
  unordered_map<int,vector<int>> particle_positions;
  {
    mt19937 random_number_generator;
    random_number_generator.seed(seed+1);
    uniform_int_distribution<int> distribution(0,distance-1);
    int particle_index = 0;
    while(particle_index<particles){
      int x = 0;
      assert(x<distance);
      assert(x>=0);
      int siteId = converter.to1D(x,y,z);
      if(siteOccupied.count(siteId)==0){
        vector<int> position = { x, y, z};
        particle_positions[particle_index] = position;
        ++particle_index;
        siteOccupied.insert(siteId);
      }
    }
  } // Place particles randomly in the system  

  // Run course grained Monte Carlo
  cout << "Running course grained Monte Carlo" << endl;
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

    class Electron : public KMC_Particle {};
    // Create the electrons using the KMC_Particle class
    vector<KMC_Particle> electrons;        
    {
      for(int particle_index = 0; particle_index<particles; ++particle_index){
        Electron electron;
        int siteId = converter.to1D(particle_positions[particle_index]);
        electron.occupySite(siteId);
        electrons.push_back(electron);
      }
    }
    
    // Run the course grain simulation
    {
      KMC_CourseGrainSystem CGsystem;
      CGsystem.setRandomSeed(seed);
      CGsystem.setMinCourseGrainIterationThreshold(threshold);
      CGsystem.setMaxCourseGrainResolution(resolution);
      CGsystem.initializeSystem(rates_to_neighbors);
      CGsystem.initializeParticles(electrons);

      double current_time_sample_increment = cutoff_time/static_cast<double>(sample_rate);
      double sample_time = current_time_sample_increment;

      vector<double> transient_current(sample_rate,0.0);
      // Calculate Particle dwell times and sort 
      list<pair<int,double>> particle_global_times;
      {

        mt19937 random_number_generator;
        random_number_generator.seed(3);
        uniform_real_distribution<double> distribution(0.0,1.0);

        for(int particle_index=0; particle_index<particles;++particle_index){
          particle_global_times.push_back(pair<int,double>(particle_index,electrons.at(particle_index).getDwellTime()));
        }
        particle_global_times.sort(compareSecondItemOfPair);
      }// Calculate particle dwell times and sort
      assert(particle_global_times.begin()->second<cutoff_time);

   
      int current_index = 0; 
      while(particle_global_times.size()>0 && particle_global_times.begin()->second<cutoff_time){
        double deltaX = 0.0;
        while(particle_global_times.size()>0 && particle_global_times.begin()->second<sample_time){
          auto particle_index = particle_global_times.begin()->first;
          KMC_Particle& electron = electrons.at(particle_index); 
          int siteId = electron.getIdOfSiteCurrentlyOccupying();
          int old_x_pos = siteId;
          CGsystem.hop(electron);
          siteId = electron.getIdOfSiteCurrentlyOccupying();
          int new_x_pos = siteId;
          deltaX+=static_cast<double>(new_x_pos-old_x_pos);
          // Update the dwell time
          particle_global_times.begin()->second += electron.getDwellTime();
          // reorder the particles based on which one will move next
          if(new_x_pos==(distance-1)){
            CGsystem.removeParticleFromSystem(electron);
            particle_global_times.pop_front();
          }
          particle_global_times.sort(compareSecondItemOfPair);
        }

        transient_current.at(current_index) = deltaX/current_time_sample_increment;
        ++current_index;
        sample_time+=current_time_sample_increment;
      }

      cout << endl;
      cout << "Transient Current" << endl;
      sample_time = 0.0;
      for( auto current : transient_current){
        sample_time+=current_time_sample_increment;
        cout << sample_time << " " << current << endl;
      }  

    }// End of the Course grain simulation 
    

  } // End of course grain Monte Carlo

  return 0;
}
