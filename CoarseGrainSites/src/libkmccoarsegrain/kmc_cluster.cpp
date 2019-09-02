
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <cassert>

#include "../../include/kmccoarsegrain/kmc_constants.hpp"
#include "kmc_cluster.hpp"
#include "log.hpp"

using namespace std;

namespace kmccoarsegrain {

/****************************************************************************
 * Constants
 ****************************************************************************/

/// Cluster Id counter is used to ensure that each new cluster has a unique id
static int clusterIdCounter = 0;

/****************************************************************************
 * Public Facing Functions
 ****************************************************************************/
void KMC_Cluster::occupy(const int& siteId){
//  Site & site = sitesInCluster_[siteId];
	
  assert(site_visits_.count(siteId));
  ++total_visit_freq_; 
//  site.occupied->insert(site.id);//setToOccupiedStatus(); 
	system_ptrs_.occupied->insert(siteId);
}

void KMC_Cluster::setRandomSeed(const unsigned long seed){
  random_engine_ = mt19937(seed);
}

bool KMC_Cluster::isOccupied() const {
  for( const int & siteId : sitesInCluster_ ){
    if(system_ptrs_.occupied->count(siteId)){
      return true;
    }
  }
  return false;
}

bool KMC_Cluster::isOccupied(const int& siteId) const {
  assert(siteIsInCluster(siteId));
  return system_ptrs_.occupied->count(siteId);
}

void KMC_Cluster::vacate(const int& siteId){
  system_ptrs_.occupied->erase(siteId);
}

void KMC_Cluster::removeWalker(const int & walker_id){
  remaining_walker_dwell_times_.erase(walker_id);
}

KMC_Cluster::KMC_Cluster(KMC_CoarseGrainSystem & CGSystem) {
  id_ = clusterIdCounter;
  clusterIdCounter++;
  convergence_method_ = converge_by_iterations_per_site;
	system_ptrs_.visit_freq = &(CGSystem.crude_.visit_freq_) ;
	system_ptrs_.time_constants = &(CGSystem.crude_.time_constants_);
	system_ptrs_.occupied = &(CGSystem.crude_.site_occupied_);
	system_ptrs_.cpd_neighbors = &(CGSystem.crude_.cpd_neighbor_hop_);
	system_ptrs_.prob_neighbors = &(CGSystem.site_neigh_prob);
	system_ptrs_.neigh_rates = (CGSystem.crude_.rates_);
}

KMC_Cluster::KMC_Cluster(const KMC_Cluster & cluster) {
  id_ = cluster.id_;
  total_visit_freq_ = cluster.total_visit_freq_;
  occupied_ = cluster.occupied_;
  escape_time_constant_ = cluster.escape_time_constant_;
  random_engine_ = cluster.random_engine_;
  random_distribution_ = cluster.random_distribution_;
  prev_total_visit_freq_ = cluster.prev_total_visit_freq_;
  resolution_ = cluster.resolution_;
  iterations_ = cluster.iterations_;
  convergenceTolerance_ = cluster.convergenceTolerance_;
  convergence_method_ = cluster.convergence_method_;
  internal_time_constant_ = cluster.internal_time_constant_;
  probabilityHopToNeighbor_ = cluster.probabilityHopToNeighbor_;
  cumulitive_probabilityHopToNeighbor_ = cluster.cumulitive_probabilityHopToNeighbor_;
  internal_dwell_time_ = cluster.internal_dwell_time_;
  site_visits_ = cluster.site_visits_;
  sumOfEscapeRateFromSiteToNeighbor_ = cluster.sumOfEscapeRateFromSiteToNeighbor_;
  sumOfEscapeRateFromSiteToInternalSite_ = cluster.sumOfEscapeRateFromSiteToInternalSite_;
  sitesInCluster_ = cluster.sitesInCluster_;
  probabilityHopOffInternalSite_ = cluster.probabilityHopOffInternalSite_;
  probabilityHopBetweenInternalSite_ = cluster.probabilityHopBetweenInternalSite_;
  probabilityOnSite_ = cluster.probabilityOnSite_;
  probabilityHopToInternalSite_ = cluster.probabilityHopToInternalSite_;
  cumulitive_probabilityHopToInternalSite_ = cluster.cumulitive_probabilityHopToInternalSite_;


	system_ptrs_.visit_freq = cluster.system_ptrs_.visit_freq;
	system_ptrs_.time_constants = cluster.system_ptrs_.time_constants;
	system_ptrs_.occupied = cluster.system_ptrs_.occupied;
	system_ptrs_.cpd_neighbors = cluster.system_ptrs_.cpd_neighbors;
	system_ptrs_.prob_neighbors = cluster.system_ptrs_.prob_neighbors;
	system_ptrs_.neigh_rates = cluster.system_ptrs_.neigh_rates;
}
void KMC_Cluster::addSite(const int& siteId) {
	
	assert(siteId>-1 && "addSite error site id is negative");
	if(siteIsInCluster(siteId)){
		return;
	}
	sitesInCluster_.push_back(siteId);
	sort(sitesInCluster_.begin(),sitesInCluster_.end());
}

bool KMC_Cluster::siteIsInCluster(const int & siteId) const {
	if(binary_search(sitesInCluster_.begin(),sitesInCluster_.end(),siteId)){
		return false;	
	}
	return true;
}

// Must use copy to use std unique
void KMC_Cluster::addSites(vector<int> siteIds) {

		// Must sort before calling unique
		sort(siteIds.begin(),siteIds.end());
		unique(siteIds.begin(),siteIds.end());

		// Cannot append sites directly because will mess up binary search
		// must evaluate all sites first before appending
		vector<int> fresh_sites;
    for (const int & siteId : siteIds) {
			if(siteIsInCluster(siteId)){
				continue;
			}
			fresh_sites.push_back(siteId);
  	}
		
    sitesInCluster_.insert(sitesInCluster_.end(),fresh_sites.begin(),fresh_sites.end());
		sort(sitesInCluster_.begin(),sitesInCluster_.end());
}

void KMC_Cluster::updateProbabilitiesAndTimeConstant() {

  unordered_map<int,int> temporary_visit_frequencies = getVisitFrequencies_();

  solveMasterEquation_();
  calculateSumOfEscapeRatesFromSitesToTheirNeighbors_();
  calculateSumOfEscapeRatesFromSitesToInternalSites_();
  calculateProbabilityHopOffInternalSite_();
  calculateProbabilityHopBetweenInternalSite_();
  calculateEscapeTimeConstant_();
  calculateInternalTimeConstant_();

  site_visits_.clear();
  for( const int & siteId : sitesInCluster_ ){
    site_visits_[siteId] = 0.0;
    if(temporary_visit_frequencies.count(siteId)){
      setVisitFrequency(temporary_visit_frequencies[siteId],siteId);
    }
  }
  calculateInternalDwellTimes_();

}

unordered_map<int,int> KMC_Cluster::getVisitFrequencies_(){
  unordered_map<int,int> frequencies;

  for(const pair<int,int> & site :site_visits_){
    frequencies[site.first] = getVisitFrequency(site.first);    
  }
  return frequencies;
}

vector<int> KMC_Cluster::getSiteIdsInCluster() const {
	return sitesInCluster_;
}

vector<int> KMC_Cluster::getSiteIdsNeighboringCluster() const {
  vector<int> neighborIds;
  for(const pair<int,double> & neigh_it : probabilityHopToNeighbor_) neighborIds.push_back(neigh_it.first);
  return neighborIds;
}

double KMC_Cluster::getProbabilityOfOccupyingInternalSite(const int siteId) {
  assert(siteIsInCluster(siteId) && "the provided site is not in the cluster");
  return probabilityOnSite_[siteId];
}

void KMC_Cluster::migrateSitesFrom(KMC_Cluster& cluster) {

  unordered_map<int,int> visits;
	for( const int & siteId : cluster.sitesInCluster_ ){
    visits[siteId] = cluster.getVisitFrequency(siteId);
  }

	addSites(cluster.sitesInCluster_);

  // Change the cluster so that it will not be used unless sites are added 
  cluster.sitesInCluster_.clear();
  cluster.probabilityOnSite_.clear();
  cluster.sumOfEscapeRateFromSiteToNeighbor_.clear();
  cluster.site_visits_.clear();
  cluster.internal_dwell_time_.clear();
  cluster.probabilityHopToNeighbor_.clear();
  cluster.cumulitive_probabilityHopToNeighbor_.clear();
  cluster.escape_time_constant_ = constants::unassigned_value;
  cluster.internal_time_constant_ = constants::unassigned_value;

  updateProbabilitiesAndTimeConstant();
  // Do not need to update the cluster because the sites have been removed
  
  // Add the visits back in
  for( const pair<int,int> & site : visits){
    setVisitFrequency(site.second,site.first);
  }
}

int KMC_Cluster::pickNewSiteId(const int & walker_id) {
  if (hopWithinCluster_(walker_id)) {
		//cout << "Hopping within cluster " << endl;
    return pickInternalSite_();
  }
		//cout << "Hopping off cluster " << endl;
  return pickClusterNeighbor_(walker_id);
}

double KMC_Cluster::getProbabilityOfHoppingToNeighborOfCluster(
    const int neighId) {
  
  auto it = find_if(probabilityHopToNeighbor_.begin(),
                    probabilityHopToNeighbor_.end(),
                    [&neighId](const pair<int,double> neigh_and_prob){
                    return neigh_and_prob.first==neighId;});
  assert(it!=probabilityHopToNeighbor_.end() &&
    "Cannot get probability of hopping to neighbor, either the"
    " site is not a neighbor or the cluster has not been converged.");

  return it->second;
}

void KMC_Cluster::setConvergenceTolerance(double tolerance) {
  assert(tolerance >= 0.0 && "tolerance must be a positive value");
  convergenceTolerance_ = tolerance;
}

void KMC_Cluster::setConvergenceIterations(long iterations) {
  assert(iterations>0 && "number of iterations must be greater than 0.");
  iterations_ = iterations;
}

double KMC_Cluster::getDwellTime(const int & walker_id) {
  assert(escape_time_constant_!=constants::unassigned_value && "Cannot get "
      "dwell time of the cluster as the escape_time_constant is not defined.");
  if(remaining_walker_dwell_times_.count(walker_id)==0){
		double number = random_distribution_(random_engine_);                       
		remaining_walker_dwell_times_[walker_id]= (-1.0)*log(number) * escape_time_constant_; 
  }
  auto dwell_time = remaining_walker_dwell_times_[walker_id];
  remaining_walker_dwell_times_[walker_id]-=time_increment_;

  if(dwell_time>time_increment_){
    return time_increment_;
  }

  assert(dwell_time>0 && "Dwell time is less than 0 this means the walker "
      "should have already been removed from the cluster");
  return dwell_time;
}

void KMC_Cluster::setVisitFrequency(int frequency,const int & siteId){
  assert(site_visits_.count(siteId) && "Be sure to call occupy site to register"
      " a visit");
  assert(escape_time_constant_!=constants::unassigned_value && "Cannot set the "
      "visit frequency as the escape_time_constant is not defined. Be sure "
      "that you have called the update function and that there exist at least "
      "one rate off the cluster.");

  // Need to convert to the right storage format 
  double visits = static_cast<double>(frequency);
  site_visits_[siteId] = visits;
}

int KMC_Cluster::getVisitFrequency(const int & siteId){
  assert(site_visits_.count(siteId));
  assert(escape_time_constant_!=constants::unassigned_value && "Cannot get the "
      "visit frequency as the escape_time_constant is not defined. Be sure "
      "that you have called the update function and that there exist at least "
      "one rate off the cluster.");

  if(total_visit_freq_!=prev_total_visit_freq_){
    double difference = total_visit_freq_ - prev_total_visit_freq_; 
    for(const pair<int,int> & site_visit : site_visits_){
      double visits = static_cast<double>(difference)*probabilityOnSite_[site_visit.first]*(system_ptrs_.time_constants->at(site_visit.first));
      visits = visits/internal_dwell_time_.at(site_visit.first);
      visits = escape_time_constant_*visits;
      visits = visits/resolution_;
      site_visits_[site_visit.first] += round(visits);
    }
    prev_total_visit_freq_ = total_visit_freq_;
  }

  double visit_count = (site_visits_[siteId]); 
  return static_cast<int>(round(visit_count)); 
}

std::ostream& operator<<(std::ostream& os,
                         const kmccoarsegrain::KMC_Cluster& cluster) {

  os << "Cluster Id: " << cluster.getId() << endl;
  os << "Cluster visitFreq: " << cluster.total_visit_freq_ << endl;
  os << "Number of sites in Cluster: " << cluster.sitesInCluster_.size();
  os << endl;

  os << "Sites in cluster: " << endl;
	for( const int & siteId : cluster.sitesInCluster_ ){
		os << siteId << " ";
	}
	os << endl;
  return os;
}

double KMC_Cluster::getFastestRateOffCluster() const {
  unordered_map<int,unordered_map<int,double>> rates_to_neighbors = 
    getRatesToNeighborsOfCluster_();
  double fastest_rate = 0.0;
  for( const pair<int,unordered_map<int,double>> & site : rates_to_neighbors){
    for( const pair<int,double> & rate : site.second ) {
      if(rate.second>fastest_rate){
        fastest_rate = rate.second;
      }
    }
  }
  return fastest_rate;
}

/****************************************************************************
 * Private Internal Functions
 ****************************************************************************/

// The firsts int is the id of site within cluster
// second int is the id of the neighbor of site outside of cluster
unordered_map<int, unordered_map<int, double>>
    KMC_Cluster::getRatesToNeighborsOfCluster_() const {

  unordered_map<int, unordered_map<int, double>> externalRates;

  for (const int & siteId : sitesInCluster_) {
    for ( pair<int,double> neigh_rate : (system_ptrs_.neigh_rates->at(siteId))) {
      if (!siteIsInCluster(neigh_rate.first)) {
        //externalRates[site.first][neigh_rate.first] = neigh_rate.second;
        (*(system_ptrs_.neigh_rates))[siteId][neigh_rate.first] = neigh_rate.second;
      }
    }
  }
  return externalRates;
}

unordered_map<int, unordered_map<int,double>>
KMC_Cluster::getRatesBetweenInternalSites_() const {

  unordered_map<int, unordered_map<int,double>> internal_rates;

  for (const int & siteId : sitesInCluster_ ){
    for( const pair<int,double> & neigh_rate : (system_ptrs_.neigh_rates->at(siteId))){
      if(siteIsInCluster(neigh_rate.first)){
        internal_rates[siteId][neigh_rate.first] = neigh_rate.second;
      }
    }
  }
  return internal_rates;
}

void KMC_Cluster::initializeProbabilityOnSites_() {
  for (const int & siteId : sitesInCluster_) {
    probabilityOnSite_[siteId] =
        1.0 / (static_cast<double>(sitesInCluster_.size()));
  }
  return;
}

void KMC_Cluster::iterate_() {

  unordered_map<int, double> temp_probabilityOnSite;
  // Initialize temporary probabilities
  for(const int & siteId : sitesInCluster_){
    temp_probabilityOnSite[siteId] = probabilityOnSite_[siteId];
  }
  
  double total = 0.0;
  for (const int & siteId : sitesInCluster_) {
    for (const pair<int,double> & neigh_rate : (system_ptrs_.neigh_rates->at(siteId))) {
      if (siteIsInCluster(neigh_rate.first)) {
        temp_probabilityOnSite[siteId] +=
					system_ptrs_.prob_neighbors->at(neigh_rate.first).at(siteId) *
          probabilityOnSite_[neigh_rate.first];
      }
      temp_probabilityOnSite[siteId] -=
				system_ptrs_.prob_neighbors->at(siteId).at(neigh_rate.first) *
        probabilityOnSite_[siteId];
    }
    total += temp_probabilityOnSite[siteId];
  }

  // Combine the former probability with the presently calculated probability
  double total2 = 0.0;
  double inverse_total = 1.0/total;
  for (const int & siteId : sitesInCluster_) {
    probabilityOnSite_[siteId] =
        (temp_probabilityOnSite[siteId]*inverse_total +
         probabilityOnSite_[siteId]) * 0.5;;

    total2 += probabilityOnSite_[siteId];
  }

  // Normalize the probability
  double inverse_total2 = 1.0/total2;
  for (const int & siteId : sitesInCluster_) {
    probabilityOnSite_[siteId] = probabilityOnSite_[siteId]*inverse_total2;
  }
}

void KMC_Cluster::solveMasterEquation_() {

  initializeProbabilityOnSites_();

  if (convergence_method_ == converge_by_iterations_per_cluster) {
    for (long i = 0; i < iterations_; i++) {
      iterate_();
    }
  } else if (convergence_method_ == converge_by_iterations_per_site) {

    long total_iterations =
        iterations_ * static_cast<long>(sitesInCluster_.size());

    for (long i = 0; i < total_iterations; i++) {
      iterate_();
    }
  } else {
    double error = convergenceTolerance_ * 1.1;

    while (error > convergenceTolerance_) {
      unordered_map<int,double> oldSiteProbs = probabilityOnSite_;
      iterate_();
      error = 0.0;
      for (const pair<int,double> & site : oldSiteProbs) {
        auto diff = oldSiteProbs[site.first] - probabilityOnSite_[site.first];
        error += pow(diff, 2.0);
      }
      error = pow(error, 1.0 / 2.0);
    }
  }
  calculateProbabilityHopToInternalSite_();
  calculateProbabilityHopToNeighbors_();
  calculateInternalDwellTimes_();
}


unordered_map<int, vector<pair<int, double>>>
    KMC_Cluster::getInternalRatesFromNeighborsComingToSite_() {

  unordered_map<int, vector<pair<int, double>>> internalRates;

  for (const int & siteId : sitesInCluster_) {
    for (const pair<int,double> & neigh_rate : system_ptrs_.neigh_rates->at(siteId)) {
      if (siteIsInCluster(neigh_rate.first)) {
        pair<int, double> rateToSite(siteId, neigh_rate.second);
        internalRates[neigh_rate.first].push_back(rateToSite);
      }
    }
  }
  return internalRates;
}

bool KMC_Cluster::hopWithinCluster_(const int & walker_id) const {
  assert(remaining_walker_dwell_times_.count(walker_id) && "Walker is not "
      "found within the cluster and does not have a dwell time, error in "
      "hopWithinCluster function call. Make sure you call getDwellTime first.");
  return remaining_walker_dwell_times_.at(walker_id)>0;
}

int KMC_Cluster::pickClusterNeighbor_(const int & walker_id) {
  remaining_walker_dwell_times_.erase(walker_id);

  double number = random_distribution_(random_engine_);
 
  for(auto it = cumulitive_probabilityHopToNeighbor_.rbegin();
      it!=cumulitive_probabilityHopToNeighbor_.rend();++it){
    if (number > it->second) {
			return it->first;
		}
  }
  assert("Cummulitive probability distribution is flawed or random "
      "number is greater than 1");
  return -1;
}

int KMC_Cluster::pickInternalSite_() {

  double number = random_distribution_(random_engine_);
  for(auto it = cumulitive_probabilityHopToInternalSite_.rbegin();
      it!=cumulitive_probabilityHopToInternalSite_.rend();++it){
    if(number>it->second){
      return it->first;
    }
  }
  assert("cumulitive probability distribution of sites in the cluster "
    "is flawed or random number is greater than 1");
  return -1;
}

void KMC_Cluster::calculateProbabilityHopOffInternalSite_() {
 
  probabilityHopOffInternalSite_.clear();
  assert(sitesInCluster_.size()>1 && "Cannot create a cluster from a single site");

  unordered_map<int,unordered_map<int,double>> ratesToNeighbors = 
    getRatesToNeighborsOfCluster_();
  auto sum_rates_off = 0.0;
  for(const pair<int,unordered_map<int,double>> & site_to_neigh : ratesToNeighbors){
    for(const pair<int,double> & rate : site_to_neigh.second ){
      sum_rates_off+=rate.second;
    }
  }
  auto sum_time_constants = 0.0;
  for(const pair<int,double> & site_prob : probabilityOnSite_) {
    sum_time_constants+=(system_ptrs_.time_constants->at(site_prob.first));
  }
  double total2 = 0.0;
  for(const pair<int,double> & site_prob : probabilityOnSite_ ){
    probabilityHopOffInternalSite_[site_prob.first] = site_prob.second*sumOfEscapeRateFromSiteToNeighbor_[site_prob.first]/sum_rates_off*(system_ptrs_.time_constants->at(site_prob.first))/sum_time_constants;
    total2+=probabilityHopOffInternalSite_[site_prob.first];

  }

  for(const pair<int,double> & hop_off : probabilityHopOffInternalSite_){
    probabilityHopOffInternalSite_[hop_off.first] = hop_off.second/total2;
  }
}


void KMC_Cluster::calculateProbabilityHopBetweenInternalSite_() {
 
  probabilityHopBetweenInternalSite_.clear();
  assert(sitesInCluster_.size()>1 && "Cannot create a cluster from a single site");

  unordered_map<int,unordered_map<int,double>> ratesBetweenSites = 
    getRatesBetweenInternalSites_();

  unordered_map<int,double> sum_sites_prob_to_hop;
  double sum_internal = 0.0;

  // rate_1 to 2 / sum( rate_1 to j) is the same as rate_1 to 2 * dwell_1
  for(const pair<int,unordered_map<int,double>> & site_to_site : ratesBetweenSites){
    int site_id = site_to_site.first;
    double sum_internal_rates = 0.0;

    for(const pair<int,double> & internal_neigh : site_to_site.second){
      sum_internal_rates+=internal_neigh.second;
    }
    sum_sites_prob_to_hop[site_id] = sum_internal_rates*(system_ptrs_.time_constants->at(site_id));
    probabilityHopBetweenInternalSite_[site_id] = sum_sites_prob_to_hop[site_id]*probabilityOnSite_[site_id];
    sum_internal+=probabilityHopBetweenInternalSite_[site_id]; 
  }

  // Normalize
  for(const pair<int,double> & site : probabilityHopBetweenInternalSite_){
    probabilityHopOffInternalSite_[site.first]/=sum_internal;
  }
  
}

// Calculates the sum of the rates off of each site that is in the cluster to
// sites external to the cluster
void KMC_Cluster::calculateSumOfEscapeRatesFromSitesToTheirNeighbors_() {

  unordered_map<int,unordered_map<int,double>> ratesToNeighbors = 
    getRatesToNeighborsOfCluster_();

  unordered_map<int, double> temp;

  for (const pair<int,unordered_map<int,double>> & site : ratesToNeighbors) {
    auto site_id_source = site.first;
    for (const pair<int,double> & neigh_terminate_rate : site.second) {
      if (temp.count(site_id_source)==0) {
        temp[site_id_source] = neigh_terminate_rate.second;
      } else {
        temp[site_id_source] += neigh_terminate_rate.second;
      }
    }
  }
  sumOfEscapeRateFromSiteToNeighbor_ = temp;
}

// Calculates the sum of the rates off of each site that is in the cluster to
// sites internal to the cluster
void KMC_Cluster::calculateSumOfEscapeRatesFromSitesToInternalSites_() {

  unordered_map<int,unordered_map<int,double>> ratesToInternalSites = 
    getRatesBetweenInternalSites_();

  unordered_map<int, double> temp;
  for (const pair<int,unordered_map<int,double>> & site : ratesToInternalSites) {
    auto site_id_source = site.first;
    double sum_rates = 0.0;
    for (const pair<int,double> & site_terminate_rate : site.second) {
      sum_rates+=site_terminate_rate.second;
    }
    temp[site_id_source] += sum_rates;
  }
  sumOfEscapeRateFromSiteToInternalSite_ = temp;
}

// Requires that calculateProbabilityHopOffInternalSites has first been called
void KMC_Cluster::calculateEscapeTimeConstant_() {
  escape_time_constant_ = 0.0;
  if(getRatesToNeighborsOfCluster_().size()==0){
    escape_time_constant_ = constants::unassigned_value;
  }else{
    for( const pair<int,double> & site_prob : probabilityHopOffInternalSite_ ){
      double rate_off = sumOfEscapeRateFromSiteToNeighbor_[site_prob.first];
      if(rate_off>0){
        escape_time_constant_ += 1.0/rate_off *site_prob.second ;
      }
    }
  }
  time_increment_ = escape_time_constant_/resolution_;
}

void KMC_Cluster::calculateInternalTimeConstant_() {
  internal_time_constant_ = 0.0;
  if(getRatesBetweenInternalSites_().size()==0){
    internal_time_constant_ = constants::unassigned_value;
  }else{
    for( const pair<int,double> & site_prob : probabilityHopBetweenInternalSite_  ){
      auto rate_off = sumOfEscapeRateFromSiteToInternalSite_[site_prob.first];
      if(rate_off>0){
        internal_time_constant_ += 1.0/rate_off *site_prob.second;
      }
    }
  }
}
void KMC_Cluster::calculateProbabilityHopToInternalSite_() {

  unordered_map<int, double> probabilityHopToInternalSite;
  double total = 0.0;
  for(const pair<int,double> & site_prob : probabilityOnSite_){
    probabilityHopToInternalSite[site_prob.first] = site_prob.second * (system_ptrs_.time_constants->at(site_prob.first));

    total+=probabilityHopToInternalSite[site_prob.first];
  }

  // Normalize
  for(const pair<int,double> & site_prob_per_time : probabilityHopToInternalSite){
    probabilityHopToInternalSite[site_prob_per_time.first]/=total;
  }

  probabilityHopToInternalSite_.clear();
  copy(probabilityHopToInternalSite.begin(),
      probabilityHopToInternalSite.end(),
      back_inserter(probabilityHopToInternalSite_));

  sort(probabilityHopToInternalSite_.begin(),
      probabilityHopToInternalSite_.end(),
      [](const pair<int,double>& x,const pair<int,double>&y)->bool{
        return x.second<y.second;
      });

  total = 0.0;
  double value = 0.0;
  for(pair<int,double> site_and_prob : probabilityHopToInternalSite_){
    site_and_prob.second+=total;
    total = site_and_prob.second;
    cumulitive_probabilityHopToInternalSite_.push_back(pair<int,double>(site_and_prob.first,value));
    value = site_and_prob.second;
  }

}

// requires master equation convergence as it uses probabilityOnSite_
void KMC_Cluster::calculateProbabilityHopToNeighbors_() {

  probabilityHopToNeighbor_.clear();
  unordered_map<int, double> temp_probabilityHopToNeighbor;

  for (const int & siteId : sitesInCluster_) {
    for (const pair<int,double> & neigh_rate : (system_ptrs_.neigh_rates->at(siteId))) {
      if (!siteIsInCluster(neigh_rate.first)) {
        if(temp_probabilityHopToNeighbor.count(neigh_rate.first)){
          temp_probabilityHopToNeighbor[neigh_rate.first] +=
            system_ptrs_.prob_neighbors->at(siteId).at(neigh_rate.first) *
            probabilityOnSite_[siteId];
        }else{
          temp_probabilityHopToNeighbor[neigh_rate.first] =
						system_ptrs_.prob_neighbors->at(siteId).at(neigh_rate.first) *
            probabilityOnSite_[siteId];
        }
      }
    }
  }

  double total = 0.0;
  for(const pair<int,double> & prob : temp_probabilityHopToNeighbor){
    total += prob.second;
  }

  // Normalize the probability
  double inverse_total = 1.0/total;
  unordered_map<int, double> probabilityHopToNeighbor;
  for (const pair<int,double> & prob : temp_probabilityHopToNeighbor) {
    probabilityHopToNeighbor[prob.first] = prob.second*inverse_total;
  }

  copy(probabilityHopToNeighbor.begin(),
      probabilityHopToNeighbor.end(),
      back_inserter(probabilityHopToNeighbor_));

  sort(probabilityHopToNeighbor_.begin(),
      probabilityHopToNeighbor_.end(),
      [](const pair<int,double>& x,const pair<int,double>&y)->bool{
        return x.second<y.second;
      });

  total = 0.0;
  double value = 0.0;
  for(pair<int,double> site_and_prob : probabilityHopToNeighbor_){
    site_and_prob.second+=total;
    total = site_and_prob.second;
    cumulitive_probabilityHopToNeighbor_.push_back(pair<int,double>(site_and_prob.first,value));
    value = site_and_prob.second;
  }
}

void KMC_Cluster::calculateInternalDwellTimes_(){
  const auto & internal_rates = getRatesBetweenInternalSites_();

  for (const auto & rate_to_neigh : internal_rates) {
    double sum_internal_site_rates = 0.0;
    for (const pair<int,double> & rate : rate_to_neigh.second) {
      sum_internal_site_rates += rate.second;
    }
    internal_dwell_time_[rate_to_neigh.first] = 1.0/sum_internal_site_rates;
  }
}


}
