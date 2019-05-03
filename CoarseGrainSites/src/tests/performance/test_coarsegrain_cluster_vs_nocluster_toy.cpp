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

int main(int argc, char* argv[]){

  cout << "Runs the following toy system" << endl;
  cout << endl;
  cout << "          site1   site2         " << endl;
  cout << "            ^       ^           " << endl;
  cout << "            |       |           " << endl; 
  cout << " site3 <- site4 - site5 -> site6" << endl;
  cout << "            |       |           " << endl; 
  cout << "            v       v           " << endl;
  cout << "          site7   site8         " << endl;
  cout << endl;
  cout << "The average time between site4 and site5 are specified with the time1" << endl;
  cout << "argument. The rate1 = 1/time1" << endl;
  cout << "Averate hopping time between all other sites are controlled with the" << endl;
  cout << "time2 argument, and are only in the direction of indicated" << endl;
  cout << "by the arrows." << endl;
  cout << "A charge is placed on site4 or site5 and then runs until" << endl;
  cout << "the charge jumps to one of the other sites (excluding 4 and 5)" << endl;

  if(argc!=6){
    cerr << "To run the program correctly you must provide the " << endl;
    cerr << "following parameters: " << endl;
    cerr << endl;
    cerr << "time1      - double value defines the average time between " << endl;
    cerr << "             sites 4 and 5." << endl;
    cerr << "time2      - double value defines the average time between " << endl;
    cerr << "             all the remaining sites" << endl;
    cerr << "walkers    - integer value defines number of walkers." << endl;
    cerr << "threshold  - integer value defines minimum threshold of " << endl;
    cerr << "             how often the simulation will try to coarse " << endl;
    cerr << "             grain." << endl;
    cerr << "time       - a double defining incrementation time of the simulation " << endl;
    cerr << "increment    how often a measurement is made." << endl;
    cerr << endl;
    cerr << "To run:" << endl;
    cerr << endl;
    cerr << "./performance_test_coarsegraining_cluster_vs_nocluster time1 time2 walkers threshold time_increment" << endl;
    cerr << endl;
    return -1;
  }

  double time1 = stod(string(argv[1]));
  double time2 = stod(string(argv[2]));
  int walkers = stoi(string(argv[3]));
  int threshold = stoi(string(argv[4]));
  double time_inc = stod(string(argv[5]));

  cout << endl;
  cout << "Parameters passed in:" << endl;
  cout << endl;
  cout << "time1:       " << time1 << endl;
  cout << "time2:       " << time2 << endl;
  cout << "walkers:     " << walkers << endl;
  cout << "threshold:   " << threshold << endl;
  cout << "time inc:    " << time_inc << endl;
  cout << endl;
  cout << endl;

  double rate1 = 1/time1;
  double rate2 = 1/time2;
  unordered_map<int,unordered_map<int,double>> rates;
  /*          site1   site2        
   *            ^       ^ 
   *            |       |           
   * site3 <- site4 - site5 -> site6
   *            |       |           
   *            v       v
   *          site7   site8        
   */
  rates[4][1]=rate2;
  rates[4][3]=rate2;
  rates[4][7]=rate2;
  rates[4][5]=rate1;
  rates[5][4]=rate1;
  rates[5][2]=rate2;
  rates[5][6]=rate2;
  rates[5][8]=rate2;

  cout << "Running Crude Monte Carlo" << endl;
  // Coarse Grained Monte Carlo with no clustering
  high_resolution_clock::time_point nocluster_time_start = high_resolution_clock::now();
  {
    KMC_CoarseGrainSystem CGsystem;
    CGsystem.setRandomSeed(3);
    CGsystem.setMinCoarseGrainIterationThreshold(constants::inf_iterations);
    CGsystem.setTimeResolution(time_inc);
    CGsystem.initializeSystem(rates);
   // Run the coarse grain simulation for as many walkers as specified
    int startingSiteId = 4;
    for(int walker_index = 0; walker_index<walkers; ++walker_index){

      class Electron : public KMC_Walker {};
      vector<pair<int,KMC_Walker>> electrons; 
      Electron elec;
      elec.occupySite(startingSiteId);
      int electronId = 0;
      electrons.push_back(pair<int,KMC_Walker>(electronId,elec));

      CGsystem.initializeWalkers(electrons);
      // Calculate Walker dwell times and sort 
      list<pair<int,double>> walker_global_times;
      walker_global_times.push_back(pair<int,double>(0,electrons.at(0).second.getDwellTime()));

      int site=4;
      while(site==4 || site==5){
        int walker_index = walker_global_times.begin()->first;
        KMC_Walker& electron = electrons.at(walker_index).second; 
        int electron_id = electrons.at(walker_index).first;
        CGsystem.hop(electron_id,electron);
        // Update the dwell time
        walker_global_times.begin()->second += electron.getDwellTime();
        site = electron.getPotentialSite();
      }

      CGsystem.removeWalkerFromSystem(electrons.at(0));
    }
    auto clusters = CGsystem.getClusters();
    if(clusters.size()!=0){
      throw runtime_error("Error a cluster was found in the nocluster part of"
          " the test");
    }
  }// End of the Coarse grain simulation without clustering 
  high_resolution_clock::time_point nocluster_time_end = high_resolution_clock::now();

  double resolution_of_cluster = 0.0;  
  double time_increment_of_cluster = 0.0;
  cout << "Running coarse grained Monte Carlo" << endl;
  high_resolution_clock::time_point cluster_time_start = high_resolution_clock::now();
  { // Run with clusters
    // Run the coarse grain simulation for as many walkers as specified
    KMC_CoarseGrainSystem CGsystem;
    CGsystem.setRandomSeed(3);
    CGsystem.setMinCoarseGrainIterationThreshold(threshold);
    CGsystem.setTimeResolution(time_inc);
    CGsystem.initializeSystem(rates);

    int startingSiteId = 4;
    for(int walker_index = 0; walker_index<walkers; ++walker_index){

      class Electron : public KMC_Walker {};
      vector<pair<int,KMC_Walker>> electrons; 
      Electron elec;
      elec.occupySite(startingSiteId);
      electrons.push_back(pair<int,KMC_Walker>(0,elec));

      CGsystem.initializeWalkers(electrons);
      // Calculate Walker dwell times and sort 
      list<pair<int,double>> walker_global_times;
      walker_global_times.push_back(pair<int,double>(0,electrons.at(0).second.getDwellTime()));

      int site=4;
      while(site==4 || site==5){
        int walker_index = walker_global_times.begin()->first;
        KMC_Walker& electron = electrons.at(walker_index).second; 
        int electron_id = electrons.at(walker_index).first;
        CGsystem.hop(electron_id,electron);
        // Update the dwell time
        walker_global_times.begin()->second += electron.getDwellTime();
        site = electron.getPotentialSite();
      }

      CGsystem.removeWalkerFromSystem(electrons.at(0));
    } 
    auto clusters = CGsystem.getClusters();
    if(clusters.size()!=1){
      cerr << "WARNING a single cluster was not detected when coarse graining is run" << endl;
    }else{
      unordered_map<int,double> clusters_resolution = CGsystem.getResolutionOfClusters();
      resolution_of_cluster = clusters_resolution.begin()->second;
      unordered_map<int,double> clusters_time_inc = CGsystem.getTimeIncrementOfClusters();
      time_increment_of_cluster = clusters_time_inc.begin()->second; 
    }
  } // End of cluster coarse grain Monte Carlo
   
     high_resolution_clock::time_point cluster_time_end = high_resolution_clock::now();

  // Run coarse grained Monte Carlo

  auto duraction_nocluster = duration_cast<milliseconds>(nocluster_time_end-nocluster_time_start).count();
  auto duraction_coarse = duration_cast<milliseconds>(cluster_time_end-cluster_time_start).count();

  cout << "Crude Monte Carlo Run Time: " << duraction_nocluster << " ms" << endl;
  cout << "Coarse Monte Carlo Run Time: " << duraction_coarse << " ms";
  if(resolution_of_cluster==0.0){
    cout << " resolution NaN";
    cout << " time increment NaN" << endl;
  }else{
    cout << " resolution " << resolution_of_cluster;
    cout << " time increment " << time_increment_of_cluster << endl;
  }
  return 0;
}
