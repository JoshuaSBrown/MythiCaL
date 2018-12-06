
#include <algorithm>
#include <chrono>
#include <set>
#include <utility>
#include <cassert>

#include "kmc_site.hpp"
#include "../../../include/kmccoarsegrain/kmc_constants.hpp"

using namespace std;

namespace kmccoarsegrain {

  class CustomComparitor
  {
    public:
      int operator()(const pair<double,int>& lhs, const pair<double,int>& rhs)
      {
        if(lhs.first==rhs.first){
          return 0;
        }
        if(lhs.second==rhs.second){
          return lhs.first<rhs.first;
        }
        return lhs.second < rhs.second;
      }
  };

/*********************************************************************
 * Public Facing Functions
 *********************************************************************/

KMC_Site::KMC_Site()
    : KMC_TopologyFeature() {

  cluster_id_ = constants::unassignedId;
}

KMC_Site::~KMC_Site() {}

void KMC_Site::setRatesToNeighbors(unordered_map<int, double*> neighRates) {
  assert(neighRates.size()!=0 && "Sites must have at least one rate to a "
    "neighbor. Cannot set rates to neighbors with empty map.");
  for (auto neighAndRate : neighRates) {
    assert(*neighAndRate.second!=0 && "One of the rates is 0.0. You cannot "
        "set a rate to a value of 0.0 as it is meaningless.");
    neighRates_[neighAndRate.first] = neighAndRate.second;
  }
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

void KMC_Site::addNeighRate(const pair<int, double*> neighRate) {

  assert(neighRates_.count(neighRate.first)==0 && "That neighbor has already been added.");
  neighRates_[neighRate.first] = neighRate.second;
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

void KMC_Site::resetNeighRate(const pair<int, double*> neighRate) {
  neighRates_[neighRate.first] = neighRate.second;
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

vector<double> KMC_Site::getRateToNeighbors() {
  vector<double> rates;
  for (auto rate : neighRates_) rates.push_back(*(rate.second));
  return rates;
}

double KMC_Site::getRateToNeighbor(const int neighSiteId) {
  assert(neighRates_.count(neighSiteId)!=0 && "Error the site Id is not a neighbor of the site ");
  return *(neighRates_[neighSiteId]);
}

double KMC_Site::getFastestRate(){
  double rate =0.0;
  for(auto neigh_rate: neighRates_) {
    if(*(neigh_rate.second)>rate) rate = *(neigh_rate.second);
  }
  return rate;
}

vector<int> KMC_Site::getNeighborSiteIds() const {
  vector<int> neighborIds;
  for (auto neighId : neighRates_) neighborIds.push_back(neighId.first);
  return neighborIds;
}

int KMC_Site::pickNewSiteId() {
  double number = random_distribution_(random_engine_);
  double threshold = 0.0;
  for (auto pval : probabilityHopToNeighbor_) {
    threshold += pval.second;
    if (number < threshold) return pval.first;
  }
  assert("Error cummulitive probability distribution is flawed or "
      " the random number generator has calculated a value greater than 1, or"
      " the site has no neighbors");
  return -1;
}

unordered_map<int,double *> KMC_Site::getNeighborsAndRates(){
  return neighRates_;
}

double KMC_Site::getProbabilityOfHoppingToNeighboringSite(
    const int neighSiteId) 
{
  assert(neighRates_.count(neighSiteId) != 0 && "Error site "
      " is not a neighbor ");

  auto it = find_if(
      probabilityHopToNeighbor_.begin(),
      probabilityHopToNeighbor_.end(),
      [&neighSiteId](const pair<int,double>& neigh_and_prob)
      { return neighSiteId==neigh_and_prob.first; }); 

  return it->second;
}

vector<pair<int, double>> KMC_Site::getProbabilitiesAndIdsOfNeighbors() const {
  return probabilityHopToNeighbor_;
}

std::ostream& operator<<(std::ostream& os,
                         const kmccoarsegrain::KMC_Site& site) {
  os << "Site Id: " << site.getId() << endl;
  os << "Cluster Id: " << site.cluster_id_ << endl;
  os << "Total Visit Frequency: " << site.total_visit_freq_ << endl;
  os << "Escape Time Constant: " << site.escape_time_constant_ << endl;
  os << "Neighbors:Rates" << endl;
  for (auto rate_ptr : site.neighRates_) {
    os << "\t" << rate_ptr.first << ":" << *(rate_ptr.second) << endl;
  }
  os << "Neighbors:Probability hop to them" << endl;
  for (auto probability : site.probabilityHopToNeighbor_) {
    os << "\t" << probability.first << ":" << probability.second << endl;
  }
  return os;
}

/*********************************************************************
 * Private Internal Functions
 *********************************************************************/
void KMC_Site::calculateProbabilityHopToNeighbors_() {
  double sumRates = getSumOfRates_();
  set<pair<int,double>,CustomComparitor> neigh_and_prob;
  for (auto rateToNeigh : neighRates_) {
    auto values = pair<int,double> (rateToNeigh.first,(*rateToNeigh.second) / sumRates);
    neigh_and_prob.insert(values);
  }

  probabilityHopToNeighbor_.clear();
  copy(neigh_and_prob.begin(),neigh_and_prob.end(),back_inserter(probabilityHopToNeighbor_));

  
}

void KMC_Site::calculateDwellTimeConstant_() {
  auto sumRates = getSumOfRates_();
  escape_time_constant_ = 1.0 / sumRates;
}

double KMC_Site::getSumOfRates_() {
  double sum = 0.0;
  for (auto rate : neighRates_) {
    sum += *(rate.second);
  }
  return sum;
}

}
