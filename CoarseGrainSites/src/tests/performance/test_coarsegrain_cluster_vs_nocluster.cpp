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

#include "../../../include/kmccoarsegrain/kmc_constants.hpp"
#include "../../../include/kmccoarsegrain/kmc_coarsegrainsystem.hpp"
#include "../../../include/kmccoarsegrain/kmc_walker.hpp"

using namespace std;
using namespace std::chrono;
using namespace kmccoarsegrain;

/**
 * \brief class for converting 1d array to 3d and vice versa
 **/
class Converter {
  public: 
    Converter(int distance) { distance_=distance;};
    int to1D(int x,int y, int z){
      assert(x<distance_);
      assert(y<distance_);
      assert(z<distance_);
      return ( z* distance_ * distance_ ) + (y *distance_) + x;
    }
    int to1D(vector<int> position){
      return ( position.at(2)* distance_ * distance_ ) + (position.at(1) *distance_) + position.at(0);
    }
    vector<int> to3D(int index) {
      vector<int> position(3,0);
      position.at(2) = index / (distance_*distance_);
      index-=(position.at(2) *distance_*distance_);
      position.at(1) = index/distance_;
      position.at(0) = index % distance_;
      return position;
    }
  private:
    int distance_;
};

bool compareSecondItemOfPair(const pair<int,double> &x, const pair<int,double> & y){
  return x.second<y.second;
}

int main(int argc, char* argv[]){

  if(argc!=8){
    cerr << "To run the program correctly you must provide the " << endl;
    cerr << "following parameters: " << endl;
    cerr << endl;
    cerr << "sigma      - defines the width of the density of states it" << endl;
    cerr << "             must be a double. " << endl;
    cerr << "distance   - integer defines the width, length and height" << endl;
    cerr << "             of the simulation box in terms of the number" << endl;
    cerr << "             of sites. " << endl;
    cerr << "seed       - integer value defines the random number seed" << endl;
    cerr << "walkers    - integer value defines number of walkers." << endl;
    cerr << "threshold  - integer value defines minimum threshold of " << endl;
    cerr << "             how often the simulation will try to coarse " << endl;
    cerr << "             grain." << endl;
    cerr << "time       - a double defining how long the simulation " << endl;
    cerr << "             will run for." << endl;
    cerr << "sample     - how often the current is sampled with respect"<< endl;
    cerr << "rate         to the cutoff time." << endl;
    cerr << endl;
    cerr << "To run:" << endl;
    cerr << endl;
    cerr << "./performance_test_coarsegraining_cluster_vs_nocluster sigma distance seed walkers threshold time sample_rate" << endl;
    cerr << endl;
    return -1;
  }

  double sigma = stod(string(argv[1]));
  int distance = stoi(string(argv[2]));
  int seed     = stoi(string(argv[3]));
  int walkers = stoi(string(argv[4]));
  int threshold = stoi(string(argv[5]));
  double cutoff_time = stod(string(argv[6]));
  int rate = stoi(string(argv[7]));

  cout << endl;
  cout << "Parameters passed in:" << endl;
  cout << endl;
  cout << "sigma:       " << sigma << endl;
  cout << "distance:    " << distance << endl;
  cout << "seed:        " << seed << endl;
  cout << "walkers:     " << walkers << endl;
  cout << "threshold:   " << threshold << endl;
  cout << "time:        " << cutoff_time << endl;
  cout << "sample rate: " << rate << endl;
  cout << endl;

  /// Create Energies and place them in a vector
  cout << "Filling sites with energies from a guassian distribution " << endl;
  cout << "centered at 0.0." << endl;
  cout << "sigma of " << sigma << endl;
  cout << endl;
  cout << endl;

  // Record setup time
  high_resolution_clock::time_point setup_time_start = high_resolution_clock::now();
  vector<double> energies;
  {
    mt19937 random_number_generator;
    random_number_generator.seed(seed);
    normal_distribution<double> distribution(0.0,sigma);

    int totalNumberSites = distance*distance*distance;
    for(int i=0;i<totalNumberSites;++i){
      energies.push_back(distribution(random_number_generator));
    }
  }

  Converter converter(distance);

  unordered_map<int,unordered_map<int,double>> rates;
  unordered_map<int,vector<int>> neighbors;
  {

    double reorganization_energy = 0.01;
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
      for(int y=0;y<distance;++y){
        for(int z=0;z<distance;++z){
          // Define neighbors
          int xlow = x;
          int xhigh = x;

          int ylow = y;
          int yhigh = y;
          
          int zlow = z;
          int zhigh = z;

          if(xlow-1>0) --xlow;
          if(xhigh+1<distance) ++xhigh;
          if(ylow-1>0) --ylow;
          if(yhigh+1<distance) ++yhigh;
          if(zlow-1>0) --zlow;
          if(zhigh+1<distance) ++zhigh;

          int siteId = converter.to1D(x,y,z); 
          for( int x2 = xlow; x2<=xhigh; ++x2){
            for( int y2 = ylow; y2<=yhigh; ++y2){
              for( int z2 = zlow; z2<=zhigh; ++z2){

                assert(x2>=0);
                assert(x2<distance);
                assert(y2>=0);
                assert(y2<distance);
                assert(z2>=0);
                assert(z2<distance);
                int neighId = converter.to1D(x2,y2,z2);
                if(siteId!=neighId){
                  neighbors[siteId].push_back(neighId);
                  double deltaE = energies.at(neighId)-energies.at(siteId);
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

  // Place walkers randomly in the system
  set<int> siteOccupied;
  unordered_map<int,vector<int>> walker_positions;
  {
    mt19937 random_number_generator;
    random_number_generator.seed(seed+1);
    uniform_int_distribution<int> distribution(0,distance-1);
    int walker_index = 0;
    while(walker_index<walkers){
      int x = distribution(random_number_generator);
      int y = distribution(random_number_generator);
      int z = distribution(random_number_generator);
      assert(x<distance);
      assert(y<distance);
      assert(z<distance);
      assert(x>=0);
      assert(y>=0);
      assert(z>=0);
      int siteId = converter.to1D(x,y,z);
      if(siteOccupied.count(siteId)==0){
        vector<int> position = { x, y, z};
        walker_positions[walker_index] = position;
        ++walker_index;
        siteOccupied.insert(siteId);
      }
    }
  } // Place walkers randomly in the system  
  high_resolution_clock::time_point setup_time_end = high_resolution_clock::now();

  cout << "Running noclustering coarse grained Monte Carlo" << endl;
  // Coarse Grained Monte Carlo with no clustering
  high_resolution_clock::time_point nocluster_time_start = high_resolution_clock::now();
  {
    // greating map with pointer to rates
   /* unordered_map< int, unordered_map< int, double *>> rates_to_neighbors;
    {
      for(auto site_rates : rates){
        for( auto neigh_rate : site_rates.second){
          rates_to_neighbors[site_rates.first][neigh_rate.first] =&(rates[site_rates.first][neigh_rate.first]);
        }
      }
    }*/

    class Electron : public KMC_Walker {};
    // Create the electrons using the KMC_Walker class
    vector<pair<int,KMC_Walker>> electrons;        
    {
      for(int walker_index = 0; walker_index<walkers; ++walker_index){
        Electron electron;
        int siteId = converter.to1D(walker_positions[walker_index]);
        electron.occupySite(siteId);
        electrons.push_back(pair<int,KMC_Walker>(walker_index,electron));
      }
    }
    
    // Run the coarse grain simulation
    {
  
      double current_time_sample_increment = cutoff_time/static_cast<double>(rate);
      double sample_time = current_time_sample_increment;

      KMC_CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(seed);
      CGsystem.setMinCoarseGrainIterationThreshold(constants::inf_iterations);
      CGsystem.setTimeResolution(sample_time);
      CGsystem.initializeSystem(rates);
      CGsystem.initializeWalkers(electrons);
      // Calculate Walker dwell times and sort 
      list<pair<int,double>> walker_global_times;
      {

        mt19937 random_number_generator;
        random_number_generator.seed(3);
        uniform_real_distribution<double> distribution(0.0,1.0);

        for(int walker_index=0; walker_index<walkers;++walker_index){
          walker_global_times.push_back(pair<int,double>(walker_index,electrons.at(walker_index).second.getDwellTime()));
        }
        walker_global_times.sort(compareSecondItemOfPair);
      }// Calculate walker dwell times and sort
      assert(walker_global_times.begin()->second<cutoff_time);
      while(walker_global_times.begin()->second<cutoff_time){
        auto walker_index = walker_global_times.begin()->first;
        KMC_Walker& electron = electrons.at(walker_index).second; 
        int electron_id = electrons.at(walker_index).first; 
        CGsystem.hop(electron_id,electron);
        // Update the dwell time
        walker_global_times.begin()->second += electron.getDwellTime();
        // reorder the walkers based on which one will move next
        walker_global_times.sort(compareSecondItemOfPair);
      }

      auto clusters = CGsystem.getClusters();
      if(clusters.size()!=0){
        throw runtime_error("Error a cluster was found in the nocluster part of"
            " the test");
      } 

    }// End of the Coarse grain simulation 
    

  } // End of nocluster coarse grain Monte Carlo
  
  high_resolution_clock::time_point nocluster_time_end = high_resolution_clock::now();

  // Run coarse grained Monte Carlo
  cout << "Running coarse grained Monte Carlo" << endl;
  high_resolution_clock::time_point coarse_time_start = high_resolution_clock::now();
  {
    /*// greating map with pointer to rates
    unordered_map< int, unordered_map< int, double *>> rates_to_neighbors;
    {
      for(auto site_rates : rates){
        for( auto neigh_rate : site_rates.second){
          rates_to_neighbors[site_rates.first][neigh_rate.first] =&(rates[site_rates.first][neigh_rate.first]);
        }
      }
    }*/

    class Electron : public KMC_Walker {};
    // Create the electrons using the KMC_Walker class
    vector<pair<int, KMC_Walker>> electrons;        
    {
      for(int walker_index = 0; walker_index<walkers; ++walker_index){
        Electron electron;
        int siteId = converter.to1D(walker_positions[walker_index]);
        electron.occupySite(siteId);
        electrons.push_back(pair<int,KMC_Walker>(walker_index,electron));
      }
    }
    
    // Run the coarse grain simulation
    {

      double current_time_sample_increment = cutoff_time/static_cast<double>(rate);
      double sample_time = current_time_sample_increment;

      KMC_CoarseGrainSystem CGsystem;
      CGsystem.setRandomSeed(seed);
      CGsystem.setMinCoarseGrainIterationThreshold(threshold);
      CGsystem.setTimeResolution(sample_time);
      CGsystem.initializeSystem(rates);
      CGsystem.initializeWalkers(electrons);
      // Calculate Walker dwell times and sort 
      list<pair<int,double>> walker_global_times;
      {

        mt19937 random_number_generator;
        random_number_generator.seed(3);
        uniform_real_distribution<double> distribution(0.0,1.0);

        for(int walker_index=0; walker_index<walkers;++walker_index){
          walker_global_times.push_back(pair<int,double>(walker_index,electrons.at(walker_index).second.getDwellTime()));
        }
        walker_global_times.sort(compareSecondItemOfPair);
      }// Calculate walker dwell times and sort
      assert(walker_global_times.begin()->second<cutoff_time);
      while(walker_global_times.begin()->second<cutoff_time){
        auto walker_index = walker_global_times.begin()->first;
        KMC_Walker& electron = electrons.at(walker_index).second; 
        int electron_id = electrons.at(walker_index).first; 
        CGsystem.hop(electron_id,electron);
        // Update the dwell time
        walker_global_times.begin()->second += electron.getDwellTime();
        // reorder the walkers based on which one will move next
        walker_global_times.sort(compareSecondItemOfPair);
      }

      auto clusters = CGsystem.getClusters();
      cout << "Total number of clusters found " << clusters.size() << endl;
    }// End of the Coarse grain simulation 
    
  } // End of coarse grain Monte Carlo
  high_resolution_clock::time_point coarse_time_end = high_resolution_clock::now();

  auto duraction_nocluster = duration_cast<seconds>(setup_time_end-setup_time_start+nocluster_time_end-nocluster_time_start).count();
  auto duraction_coarse = duration_cast<seconds>(setup_time_end-setup_time_start+coarse_time_end-coarse_time_start).count();

  cout << "Coarse nocluster Monte Carlo Run Time: " << duraction_nocluster << " s " << endl;
  cout << "Coarse Monte Carlo Run Time: " << duraction_coarse << " s " << endl;
  return 0;
}
