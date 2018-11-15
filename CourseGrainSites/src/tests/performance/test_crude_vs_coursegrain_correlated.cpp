#include <iostream>
#include <vector>
#include <memory>
#include <chrono>
#include <string>
#include <random>
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
    cerr << "sigma      - defins the width of the density of states it" << endl;
    cerr << "             must be a double. " << endl;
    cerr << "distance   - integer defines the width, length and height" << endl;
    cerr << "             of the simulation box in terms of the number" << endl;
    cerr << "             of sites. " << endl;
    cerr << "seed       - integer value defines the random number " << endl;
    cerr << "             generator seed." << endl;
    cerr << "resolution - integer value defines how course the couse " << endl;
    cerr << "             the approxiation will be." << endl;
    cerr << "particles  - integer value defines number of particles." << endl;
    cerr << "threshold  - integer value defines minimum threshold of " << endl;
    cerr << "             how often the simulation will try to course " << endl;
    cerr << "             grain." << endl;
    cerr << "time       - a double defining how long the simulation " << endl;
    cerr << "             will run for." << endl;
    cerr << endl;
    cerr << "To run:" << endl;
    cerr << endl;
    cerr << "./performance_test_crude_vs_coursegraining sigma distance threshold resolution" << endl;
    cerr << endl;
    return -1;
  }

  double sigma = stod(string(argv[1]));
  int distance = stoi(string(argv[2]));
  int seed = stoi(string(argv[3]));
  int resolution = stoi(string(argv[4])); 
  int particles = stoi(string(argv[5]));
  int threshold = stoi(string(argv[6]));
  double cutoff_time = stod(string(argv[7]));

  cout << endl;
  cout << "Parameters passed in:" << endl;
  cout << endl;
  cout << "sigma:      " << sigma << endl;
  cout << "distance:   " << distance << endl;
  cout << "seed:       " << seed << endl;
  cout << "resolution: " << resolution << endl;
  cout << "particles:  " << particles << endl;
  cout << "threshold:  " << threshold << endl;
  cout << "time:       " << cutoff_time << endl;
  cout << endl;

  /// Create Energies and place them in a vector
  double time = 1.0;
  cout << "Simulating time up to " << time << " seconds " << endl;
  cout << "Filling sites with energies from a guassian distribution " << endl;
  cout << "centered at 0.0." << endl;
  cout << "sigma of " << sigma << endl;
  cout << endl;


  double simulation_cutoff_time = 1.0E-4;   
  cout << "Simulation cutoff time " << simulation_cutoff_time << " seconds";
  cout << endl;
  cout << endl;

  double percentage_coorilation_seeds = 0.001;
  double correlation_radius = 3.0;
  double max_correlation = correlation_radius*1.5;
  int inverse_percentage = static_cast<int>(1.0/percentage_coorilation_seeds);
  int num_seeds = distance*distance*distance/inverse_percentage;
  cout << "Applying corrilation" << endl;
  cout << "percent of seeds " << percentage_coorilation_seeds << endl;
  cout << "total number of seeds " << num_seeds << endl;
  cout << "correlation radius " << correlation_radius << endl;
  cout << "max reach of correlation " << max_correlation << endl;
  cout << endl;

  // Record setup time
  high_resolution_clock::time_point setup_time_start = high_resolution_clock::now();

  // Assign energies from gaussian distribution
  vector<double> energies;
  {
    mt19937 random_number_generator;
    random_number_generator.seed(1);
    normal_distribution<double> distribution(0.0,sigma);

    int totalNumberSites = distance*distance*distance;
    for(int i=0;i<totalNumberSites;++i){
      energies.push_back(distribution(random_number_generator));
    }
  }

  Converter converter(distance);
  // Correlate energies
  {
    mt19937 random_number_generator;
    random_number_generator.seed(7);
    uniform_int_distribution<int> distribution(0,distance*distance*distance-1);

    vector<int> seeds;
    for(int seed_index = 0 ;seed_index< num_seeds;++seed_index){
      seeds.push_back(distribution(random_number_generator));
    }

    double distance_double = static_cast<double>(distance);
    for( auto seed : seeds ){
      energies.at(seed) = 0.0;
      auto position = converter.to3D(seed);
      double x_pos = static_cast<double>(position.at(0));
      double y_pos = static_cast<double>(position.at(1));
      double z_pos = static_cast<double>(position.at(2));

      for(double x_neigh = x_pos - max_correlation; x_neigh< x_pos+max_correlation; x_neigh+=1.0){
        for(double y_neigh = y_pos - max_correlation; y_neigh< y_pos+max_correlation; y_neigh+=1.0){
          for(double z_neigh = z_pos - max_correlation; z_neigh< z_pos+max_correlation; z_neigh+=1.0){
            if(x_neigh<distance_double && x_neigh>-0.0){
              if(y_neigh<distance_double && y_neigh>-0.0){
                if(z_neigh<distance_double && z_neigh>-0.0){
                  auto xdiff = pow(x_pos-x_neigh,2.0);
                  auto ydiff = pow(y_pos-y_neigh,2.0);
                  auto zdiff = pow(z_pos-z_neigh,2.0);
                  double neigh_radius = pow(xdiff+ydiff+zdiff,1.0/2.0);
                  int x_n = static_cast<int>(x_neigh);  
                  int y_n = static_cast<int>(y_neigh);  
                  int z_n = static_cast<int>(z_neigh);  
                  int neighId = converter.to1D(x_n,y_n,z_n);
                  assert(neighId<distance*distance*distance);
                  double energyDiff = energies.at(seed)-energies.at(neighId);
                  energies.at(neighId) = energies.at(neighId)+energyDiff*exp(-neigh_radius/correlation_radius);
                }
              } 
            }
          }
        }
      }
    }
  }

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

  // Place particles randomly in the system
  set<int> siteOccupied;
  unordered_map<int,vector<int>> particle_positions;
  {
    mt19937 random_number_generator;
    random_number_generator.seed(2);
    uniform_int_distribution<int> distribution(0,distance-1);
    int particle_index = 0;
    while(particle_index<particles){
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
        particle_positions[particle_index] = position;
        ++particle_index;
        siteOccupied.insert(siteId);
      }
    }
  } // Place particles randomly in the system  
  high_resolution_clock::time_point setup_time_end = high_resolution_clock::now();

  cout << "Running crude Monte Carlo" << endl;
  // Crude Monte Carlo
  high_resolution_clock::time_point crude_time_start = high_resolution_clock::now();
  {

    unordered_map<int,double> sojourn_times;
    unordered_map<int,double> sum_rates;
    // Calculate sojourn times & sum_rates
    {
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

            double sum_times = 0.0;
            double sum_rate = 0.0;
            int siteId = converter.to1D(x,y,z); 
            for( int x2 = xlow; x2<=xhigh; ++x2){
              for( int y2 = ylow; y2<=yhigh; ++y2){
                for( int z2 = zlow; z2<=zhigh; ++z2){
                  int neighId = converter.to1D(x2,y2,z2);
                  if(siteId!=neighId){
                    sum_times+=1.0/rates[siteId][neighId];
                    sum_rate +=rates[siteId][neighId];
                  }
                }
              }
            }
            sojourn_times[siteId] = sum_times;        
            sum_rates[siteId] = sum_rate;  
          }
        }
      }
    }// Calculate sojourn times & sum_rates

    unordered_map<int,unordered_map<int,double>> cummulitive_probability_to_neighbors;
    // Calculate crude probability to neighbors
    {
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

            unordered_map<int,double> cummulitive_probability;
            double pval = 0.0;
            int siteId = converter.to1D(x,y,z); 
            for( int x2 = xlow; x2<=xhigh; ++x2){
              for( int y2 = ylow; y2<=yhigh; ++y2){
                for( int z2 = zlow; z2<=zhigh; ++z2){
                  int neighId = converter.to1D(x2,y2,z2);
                  if(siteId!=neighId){
                    cummulitive_probability[neighId]=rates[siteId][neighId]/sum_rates[siteId];
                    cummulitive_probability[neighId]+=pval;
                    pval+=rates[siteId][neighId]/sum_rates[siteId];
                  }
                }
              }
            }
            assert(cummulitive_probability.size()!=0);
            cummulitive_probability_to_neighbors[siteId] = cummulitive_probability; 
            assert(cummulitive_probability_to_neighbors[siteId].size()!=0);
          }
        }
      }
    }// Calculate crude probability to neighbors


    // Calculate Particle dwell times and sort 
    list<pair<int,double>> particle_global_times;
    {

      mt19937 random_number_generator;
      random_number_generator.seed(3);
      uniform_real_distribution<double> distribution(0.0,1.0);

      for(int particle_index=0; particle_index<particles;++particle_index){
        auto position = particle_positions[particle_index];
        auto siteId = converter.to1D(position);
        particle_global_times.push_back(pair<int,double>(particle_index, sojourn_times[siteId]*log(distribution(random_number_generator))*-1.0));
      }
      particle_global_times.sort(compareSecondItemOfPair);
    }// Calculate particle dwell times and sort


    // Run simulation until cutoff simulation time is reached
    {

      mt19937 random_number_generator;
      random_number_generator.seed(seed);
      uniform_real_distribution<double> distribution(0.0,1.0);

      while(particle_global_times.begin()->second<cutoff_time){
        int particleId = particle_global_times.begin()->first;
        vector<int> particle_position = particle_positions[particleId];
        int siteId = converter.to1D(particle_position);

        double random_number = distribution(random_number_generator);
        // Attempt to hop
        assert(cummulitive_probability_to_neighbors[siteId].size()!=0);
        for( auto pval_iterator : cummulitive_probability_to_neighbors[siteId] ){
          if(random_number < pval_iterator.second){
            int neighId = pval_iterator.first;
            if(siteOccupied.count(neighId)){
              // Update the sojourn time particle is unable to make the jump
              particle_global_times.begin()->second += sojourn_times[siteId]*log(distribution(random_number_generator))*-1.0;
            }else{
              // vacate site
              siteOccupied.erase(siteId); 
              // Occupy new site
              siteOccupied.insert(neighId);
              // Update the particles position
              particle_positions[particleId] = converter.to3D(neighId);
              // Update the sojourn time of the particle
              particle_global_times.begin()->second += sojourn_times[neighId]*log(distribution(random_number_generator))*-1.0;
            }
            break;
          }
        }
        // reorder the particles based on which one will move next
        particle_global_times.sort(compareSecondItemOfPair);
      }
  
    } // Run simulation until cutoff simulation time is reached

    
  } // Crude Mone Carlo
  high_resolution_clock::time_point crude_time_end = high_resolution_clock::now();

  // Run course grained Monte Carlo
  cout << "Running course grained Monte Carlo" << endl;
  high_resolution_clock::time_point course_time_start = high_resolution_clock::now();
  {
    // greating map with pointer to rates
    unordered_map< int, unordered_map< int, double *>> rates_to_neighbors;
    {
      for(auto site_rates : rates){
        unordered_map< int ,double *> rates_to;
        for( auto neigh_rate : site_rates.second){
          rates_to_neighbors[site_rates.first][neigh_rate.first] = &neigh_rate.second;
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

      while(particle_global_times.begin()->second<cutoff_time){
        auto particle_index = particle_global_times.begin()->first;
        KMC_Particle & electron = electrons.at(particle_index); 
        CGsystem.hop(electron);
        // Update the dwell time
        particle_global_times.begin()->second += electron.getDwellTime();
        // reorder the particles based on which one will move next
        particle_global_times.sort(compareSecondItemOfPair);
      }

    }// End of the Course grain simulation 
    

  } // End of course grain Monte Carlo
  high_resolution_clock::time_point course_time_end = high_resolution_clock::now();

  auto duraction_crude = duration_cast<seconds>(setup_time_end-setup_time_start+crude_time_end-crude_time_start).count();
  auto duraction_course = duration_cast<seconds>(setup_time_end-setup_time_start+course_time_end-course_time_start).count();

  cout << "Crude Monte Carlo Run Time: " << duraction_crude << " s " << endl;
  cout << "Course Monte Carlo Run Time: " << duraction_course << " s " << endl;
  return 0;
}
