
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <cassert>

#include "kmc_cluster.hpp"
#include "kmc_site.hpp"
#include "../log.hpp"

using namespace std;

namespace kmccoursegrain {

/****************************************************************************
 * Constants
 ****************************************************************************/

/// Cluster Id counter is used to ensure that each new cluster has a unique id
static int clusterIdCounter = 0;

/****************************************************************************
 * Public Facing Functions
 ****************************************************************************/
void occupyCluster_(KMC_TopologyFeature* feature, int& siteId){
  auto cluster = static_cast<KMC_Cluster *>(feature);
  
  KMC_Site & site = cluster->sitesInCluster_[siteId];
  int site_visits = site.getVisitFrequency();
  site_visits += cluster->resolution_;
  site.setVisitFrequency(site_visits);
  site.setToOccupiedStatus(); 
}

bool isOccupiedCluster_(KMC_TopologyFeature* feature,int& siteId){
  auto cluster = static_cast<KMC_Cluster *>(feature);
  assert(cluster->sitesInCluster_.count(siteId));
  return cluster->sitesInCluster_[siteId].isOccupied();
}

void vacateCluster_(KMC_TopologyFeature* feature,int& siteId){
  auto cluster = static_cast<KMC_Cluster *>(feature);
  cluster->sitesInCluster_[siteId].vacate();
  cluster->vacate();
}

KMC_Cluster::KMC_Cluster() : KMC_TopologyFeature() {
  setId(clusterIdCounter);
  clusterIdCounter++;
  iterations_ = 3;
  resolution_ = 20;
  total_visit_freq_ = 0;
  convergenceTolerance_ = 0.01;
  convergence_method_ = converge_by_iterations_per_site;

  occupy_siteId_ptr_ = occupyCluster_;
  vacate_siteId_ptr_ = vacateCluster_;
  isOccupied_siteId_ptr_ = isOccupiedCluster_;

}

void KMC_Cluster::addSite(KMC_Site& newSite) {
//  if (sitesInCluster_.count(newSite.getId())) {
//    throw invalid_argument("Site has already been added to the cluster");
//  }
  assert(sitesInCluster_.count(newSite.getId())==0 && "Site has already been added to the cluster");
  newSite.setClusterId(this->getId());
  sitesInCluster_[newSite.getId()] = newSite;

  updateProbabilitiesAndTimeConstant();
}

void KMC_Cluster::addSites(vector<KMC_Site>& newSites) {
  for (KMC_Site & site : newSites) {
//    if (sitesInCluster_.count(site.getId())) {
//      throw invalid_argument("Site has already been added to the cluster");
//    }
    assert(sitesInCluster_.count(site.getId())==0 && "Site has already been added to the cluster");
    site.setClusterId(this->getId());
    sitesInCluster_[site.getId()] = site;
  }
  updateProbabilitiesAndTimeConstant();
}

void KMC_Cluster::updateProbabilitiesAndTimeConstant() {
  solveMasterEquation_();
  calculateEscapeRatesFromSitesToTheirNeighbors_();
  calculateProbabilityHopOffInternalSite_();
  calculateEscapeTimeConstant_();
}

vector<KMC_Site> KMC_Cluster::getSitesInCluster() const {
  vector<KMC_Site> sites;
  for (auto site : sitesInCluster_) sites.push_back(site.second);
  return sites;
}

vector<int> KMC_Cluster::getSiteIdsInCluster() const {
  vector<int> siteIds;
  for (auto site : sitesInCluster_) {
    siteIds.push_back(site.first);
  }
  return siteIds;
}

vector<int> KMC_Cluster::getSiteIdsNeighboringCluster() const {
  vector<int> neighborIds;
  for(auto neigh_it : probabilityHopToNeighbor_) neighborIds.push_back(neigh_it.first);
  return neighborIds;
}

double KMC_Cluster::getProbabilityOfOccupyingInternalSite(const int siteId) {
//  if (!sitesInCluster_.count(siteId)) {
//    throw invalid_argument("the provided site is not in the cluster");
//  }
  assert(sitesInCluster_.count(siteId) && "the provided site is not in the cluster");
  return probabilityOnSite_[siteId];
}

void KMC_Cluster::migrateSitesFrom(KMC_Cluster& cluster) {

  cout << "Migrating sites from one cluster to the other" << endl;
  for (auto& site : cluster.sitesInCluster_) {
    site.second.setClusterId(getId());
  }

  move(cluster.sitesInCluster_.begin(),
      cluster.sitesInCluster_.end(),
      inserter(this->sitesInCluster_,this->sitesInCluster_.end()));

  cluster.sitesInCluster_.clear();

  updateProbabilitiesAndTimeConstant();
  cluster.updateProbabilitiesAndTimeConstant();
}

int KMC_Cluster::pickNewSiteId() {
  if (hopWithinCluster_()) {
    return pickInternalSite_();
  }
  return pickClusterNeighbor_();
}

double KMC_Cluster::getProbabilityOfHoppingToNeighborOfCluster(
    const int neighId) {
  
  auto it = find_if(probabilityHopToNeighbor_.begin(),
                    probabilityHopToNeighbor_.end(),
                    [&neighId](const pair<int,double> neigh_and_prob){
                    return neigh_and_prob.first==neighId;});
//  if (it==probabilityHopToNeighbor_.end()) {
//    string err =
//        "Cannot get probability of hopping to neighbor, either the"
//        " site is not a neighbor or the cluster has not been converged.";
//    throw invalid_argument(err);
//  }
  assert(it!=probabilityHopToNeighbor_.end() &&
    "Cannot get probability of hopping to neighbor, either the"
    " site is not a neighbor or the cluster has not been converged.");

  return it->second;
}

void KMC_Cluster::setConvergenceTolerance(double tolerance) {
//  if (tolerance < 0.0) {
//    throw invalid_argument("tolerance must be a positive value");
//  }
  assert(tolerance >= 0.0 && "tolerance must be a positive value");
  convergenceTolerance_ = tolerance;
}

void KMC_Cluster::setConvergenceIterations(long iterations) {
//  if (iterations < 1) {
//    throw invalid_argument("number of iterations must be greater than 0.");
//  }
  assert(iterations>0 && "number of iterations must be greater than 0.");
  iterations_ = iterations;
}

double KMC_Cluster::getDwellTime() {
  return (KMC_TopologyFeature::getDwellTime() / resolution_);
}

void KMC_Cluster::setVisitFrequency(int frequency, int siteId){
  sitesInCluster_[siteId].setVisitFrequency(frequency);
}

int KMC_Cluster::getVisitFrequency(int siteId){
  return sitesInCluster_[siteId].getVisitFrequency();
}

std::ostream& operator<<(std::ostream& os,
                         const kmccoursegrain::KMC_Cluster& cluster) {

  os << "Cluster Id: " << cluster.getId() << endl;
  os << "Cluster visitFreq: " << cluster.total_visit_freq_ << endl;
  os << "Number of sites in Cluster: " << cluster.sitesInCluster_.size();
  os << endl;

  os << "Sites in cluster: " << endl;
  for (auto site : cluster.sitesInCluster_) {
    os << (site.second) << endl;
  }
  return os;
}

double KMC_Cluster::getFastestRateOffCluster(){
  auto rates_to_neighbors = getRatesToNeighborsOfCluster_();
  double fastest_rate = 0.0;
  for( auto site_pr : rates_to_neighbors){
    for( auto rate : site_pr.second ) {
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
    KMC_Cluster::getRatesToNeighborsOfCluster_() {

  unordered_map<int, unordered_map<int, double>> externalRates;

  for (auto site : sitesInCluster_) {
    for (auto neighId : site.second.getNeighborSiteIds()) {
      if (!siteIsInCluster(neighId)) {
        externalRates[site.first][neighId] =
            site.second.getRateToNeighbor(neighId);
      }
    }
  }
  return externalRates;
}

void KMC_Cluster::initializeProbabilityOnSites_() {
  for (auto site : sitesInCluster_) {
    probabilityOnSite_[site.first] =
        1.0 / (static_cast<double>(sitesInCluster_.size()));
  }
  return;
}

void KMC_Cluster::calculateProbabilityHopToInternalSite_() {

  unordered_map<int, double> probabilityHopToInternalSite;
  double total = 0.0;
  for(auto site_prob : probabilityOnSite_){
    probabilityHopToInternalSite[site_prob.first] = site_prob.second * sitesInCluster_[site_prob.first].getTimeConstant();

    total+=probabilityHopToInternalSite[site_prob.first];
  }

  // Normalize
  for(auto site_prob_per_time : probabilityHopToInternalSite){
    probabilityHopToInternalSite[site_prob_per_time.first]/=total;
  }

  probabilityHopToInternalSite_.clear();
  copy(probabilityHopToInternalSite.begin(),
      probabilityHopToInternalSite.end(),
      back_inserter(probabilityHopToInternalSite_));

  sort(probabilityHopToInternalSite_.begin(),
      probabilityHopToInternalSite_.end(),
      [](const pair<int,double>& x,const pair<int,double>&y)->bool{
        return x.second>y.second;
      });

}

void KMC_Cluster::calculateProbabilityHopToNeighbors_() {

  auto ratesToNeighbors = getRatesToNeighborsOfCluster_();

  double sumRatesOffCluster = 0.0;
  for (auto rateToNeigh : ratesToNeighbors) {
    for (auto rate : rateToNeigh.second) {
      sumRatesOffCluster += rate.second;
    }
  }
  auto sumDwell = 0.0;
  for (auto site : sitesInCluster_) {
    sumDwell += site.second.getTimeConstant();
  }

  unordered_map<int, double> probabilityHopToNeighbor;
  for (auto rateToNeigh : ratesToNeighbors) {
    int siteHoppingFrom = rateToNeigh.first;
    for (auto rate : rateToNeigh.second) {
      int siteHoppingTo = rate.first;

      if (probabilityHopToNeighbor.count(siteHoppingTo)) {
        probabilityHopToNeighbor[siteHoppingTo] +=
            probabilityOnSite_[siteHoppingFrom] *
            sitesInCluster_[siteHoppingFrom].getTimeConstant() / sumDwell *
            sitesInCluster_[siteHoppingFrom].getRateToNeighbor(siteHoppingTo) /
            sumRatesOffCluster;

      } else {

        probabilityHopToNeighbor[siteHoppingTo] =
            probabilityOnSite_[siteHoppingFrom] *
            sitesInCluster_[siteHoppingFrom].getTimeConstant() / sumDwell *
            sitesInCluster_[siteHoppingFrom].getRateToNeighbor(siteHoppingTo) /
            sumRatesOffCluster;
      }
    }
  }

  double total = 0.0;
  for (auto neighborProb : probabilityHopToNeighbor) {
    total += neighborProb.second;
  }
  for (auto neighborProb : probabilityHopToNeighbor) {
    probabilityHopToNeighbor[neighborProb.first] = neighborProb.second / total;
  }
  probabilityHopToNeighbor_.clear();

  copy(probabilityHopToNeighbor.begin(),
      probabilityHopToNeighbor.end(),
      back_inserter(probabilityHopToNeighbor_));

  sort(probabilityHopToNeighbor_.begin(),
      probabilityHopToNeighbor_.end(),
      [](const pair<int,double>& x,const pair<int,double>&y)->bool{
        return x.second>y.second;
      });
}

void KMC_Cluster::iterate_() {

  unordered_map<int, double> temp_probabilityOnSite;

  double total = 0.0;
  for (auto site : sitesInCluster_) {
    for (auto neighsite : site.second.getNeighborSiteIds()) {
      if (siteIsInCluster(neighsite)) {
        if (temp_probabilityOnSite.count(site.first)) {
          temp_probabilityOnSite[site.first] +=
              sitesInCluster_[neighsite].getProbabilityOfHoppingToNeighboringSite(site.first) *
              probabilityOnSite_[neighsite];
                  
        } else {
          temp_probabilityOnSite[site.first] =
              sitesInCluster_[neighsite].getProbabilityOfHoppingToNeighboringSite(site.first) *
              probabilityOnSite_[neighsite];
        }
        total += temp_probabilityOnSite[site.first];
      }
    }
  }

  double total2 = 0.0;
  for (auto site : sitesInCluster_) {
    probabilityOnSite_[site.first] =
        (temp_probabilityOnSite[site.first] / total +
         probabilityOnSite_[site.first]) /
        2.0;

    total2 += probabilityOnSite_[site.first];
  }

  for (auto site : sitesInCluster_) {
    probabilityOnSite_[site.first] = probabilityOnSite_[site.first] / total2;
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
      auto oldSiteProbs = probabilityOnSite_;
      iterate_();
      error = 0.0;
      for (auto site : oldSiteProbs) {
        auto diff = oldSiteProbs[site.first] - probabilityOnSite_[site.first];
        error += pow(diff, 2.0);
      }
      error = pow(error, 1.0 / 2.0);
    }
  }
  calculateProbabilityHopToInternalSite_();
  calculateProbabilityHopToNeighbors_();
}

unordered_map<int, vector<pair<int, double>>>
    KMC_Cluster::getInternalRatesFromNeighborsComingToSite_() {

  unordered_map<int, vector<pair<int, double>>> internalRates;

  for (auto site : sitesInCluster_) {
    for (auto neighId : site.second.getNeighborSiteIds()) {
      if (siteIsInCluster(neighId)) {
        auto rateToNeigh = site.second.getRateToNeighbor(neighId);
        pair<int, double> rateToSite(site.first, rateToNeigh);
        internalRates[neighId].push_back(rateToSite);
      }
    }
  }
  return internalRates;
}

bool KMC_Cluster::hopWithinCluster_() {
  double hop = random_distribution_(random_engine_);
  double in_or_out_threshold = (static_cast<double>(resolution_) - 1.0) /
                            static_cast<double>(resolution_);
  return (hop < in_or_out_threshold);
}

int KMC_Cluster::pickClusterNeighbor_() {

  double number = random_distribution_(random_engine_);
  double threshold = 0.0;
  for (auto pval : probabilityHopToNeighbor_) {
    threshold += pval.second;
    if (number < threshold) return pval.first;
  }
//  string err =
//      "Cummulitive probability distribution is flawed or random "
//      "number is greater than 1";
//  throw invalid_argument(err);
  assert("Cummulitive probability distribution is flawed or random "
      "number is greater than 1");
  return -1;
}

int KMC_Cluster::pickInternalSite_() {

  double number = random_distribution_(random_engine_);
  double threshold = 0.0;
  for (auto pval : probabilityHopToInternalSite_) {
    threshold += pval.second;
    if (number < threshold) {
      return pval.first;
    }
  }
  //string err =
  //    "cummulitive probability distribution of sites in the cluster "
  //    "is flawed or random number is greater than 1";
  //throw invalid_argument(err);
  assert("cummulitive probability distribution of sites in the cluster "
    "is flawed or random number is greater than 1");
  return -1;
}

void KMC_Cluster::calculateProbabilityHopOffInternalSite_() {
 
  probabilityHopOffInternalSite_.clear();

  auto ratesToNeighbors = getRatesToNeighborsOfCluster_();
  auto sum_rates_off = 0.0;
  for(auto site_to_neigh : ratesToNeighbors){
    for(auto rate : site_to_neigh.second ){
      sum_rates_off+=rate.second;
    }
  }
  auto sum_time_constants = 0.0;
  for(auto site_prob : probabilityOnSite_) {
    sum_time_constants+=sitesInCluster_[site_prob.first].getTimeConstant();
  }
  double total2 = 0.0;
  for(auto site_prob : probabilityOnSite_ ){
    probabilityHopOffInternalSite_[site_prob.first] = site_prob.second*escapeRateFromSiteToNeighbor_[site_prob.first]/sum_rates_off*sitesInCluster_[site_prob.first].getTimeConstant()/sum_time_constants;
    total2+=probabilityHopOffInternalSite_[site_prob.first];

  }

  for(auto hop_off : probabilityHopOffInternalSite_){
    probabilityHopOffInternalSite_[hop_off.first] = hop_off.second/total2;
  }
  
}

void KMC_Cluster::calculateEscapeRatesFromSitesToTheirNeighbors_() {

  auto ratesToNeighbors = getRatesToNeighborsOfCluster_();

  unordered_map<int, double> temp;

  for (auto site : ratesToNeighbors) {
    for (auto neigh : site.second) {
      if (temp.count(site.first)==0) {
        temp[site.first] = neigh.second;
      } else {
        temp[site.first] += neigh.second;
      }
    }
  }
  escapeRateFromSiteToNeighbor_ = temp;
}

void KMC_Cluster::calculateEscapeTimeConstant_() {
  auto ratesToNeighbors = getRatesToNeighborsOfCluster_();
  escape_time_constant_ = 0.0;

  for( auto site_prob : probabilityHopOffInternalSite_ ){
    escape_time_constant_ += 1.0/escapeRateFromSiteToNeighbor_[site_prob.first] *site_prob.second ;
  }
}
}
