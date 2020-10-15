
#include <algorithm>
#include <chrono>
#include <set>
#include <utility>
#include <cassert>

#include "site.hpp"
#include "mythical/constants.hpp"

using namespace std;

namespace mythical {

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

Site::Site()
    : TopologyFeature() {

  cluster_id_ = constants::unassignedId;
}

Site::~Site() {}

void Site::setRatesToNeighbors(unordered_map<int, double>& neighRates) {
  assert(neighRates.size()!=0 && "Sites must have at least one rate to a "
    "neighbor. Cannot set rates to neighbors with empty map.");
  for (auto neighAndRate : neighRates) {
    assert(neighAndRate.second!=0 && "One of the rates is 0.0. You cannot "
        "set a rate to a value of 0.0 as it is meaningless.");
    neighRates_[neighAndRate.first] = &(neighRates[neighAndRate.first]);
  }
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

void Site::addNeighRate(const pair<int, double*> neighRate) {

  assert(neighRates_.count(neighRate.first)==0 && "That neighbor has already been added.");
  neighRates_[neighRate.first] = neighRate.second;
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

void Site::resetNeighRate(const pair<int, double*> neighRate) {
  neighRates_[neighRate.first] = neighRate.second;
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

vector<double> Site::getRateToNeighbors() const {
  vector<double> rates;
  for (auto & rate : neighRates_) rates.push_back(*(rate.second));
  return rates;
}

double Site::getRateToNeighbor(const int & neighSiteId) const {
  assert(neighRates_.count(neighSiteId)!=0 && "Error the site Id is not a neighbor of the site ");
  return *(neighRates_.at(neighSiteId));
}

double Site::getFastestRate(){
  double rate =0.0;
  for(auto neigh_rate: neighRates_) {
    if(*(neigh_rate.second)>rate) rate = *(neigh_rate.second);
  }
  return rate;
}

vector<int> Site::getNeighborSiteIds() const {
  vector<int> neighborIds;
  for (auto neighId : neighRates_) neighborIds.push_back(neighId.first);
  return neighborIds;
}

int Site::pickNewSiteId(const int &) {
  return pickNewSiteId();
}

int Site::pickNewSiteId() {
  double number = random_distribution_(random_engine_);
  double threshold = 0.0;
  for (pair<int,double> & pval : probabilityHopToNeighbor_) {
    threshold += pval.second;
    if (number < threshold) return pval.first;
  }
  assert("Error cummulitive probability distribution is flawed or "
      " the random number generator has calculated a value greater than 1, or"
      " the site has no neighbors");
  return -1;
}

unordered_map<int,double *> Site::getNeighborsAndRates(){
  return neighRates_;
}

const unordered_map<int,double *> & Site::getNeighborsAndRatesConst() const{
  return neighRates_;
}

double Site::getProbabilityOfHoppingToNeighboringSite(
    const int & neighSiteId) 
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

vector<pair<int, double>> Site::getProbabilitiesAndIdsOfNeighbors() const {
  return probabilityHopToNeighbor_;
}

std::ostream& operator<<(std::ostream& os,
                         const Site& site) {
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
void Site::calculateProbabilityHopToNeighbors_() {
  double sumRates = getSumOfRates_();
  set<pair<int,double>,CustomComparitor> neigh_and_prob;
  for (auto rateToNeigh : neighRates_) {
    auto values = pair<int,double> (rateToNeigh.first,(*rateToNeigh.second) / sumRates);
    neigh_and_prob.insert(values);
  }

  probabilityHopToNeighbor_.clear();
  copy(neigh_and_prob.begin(),neigh_and_prob.end(),back_inserter(probabilityHopToNeighbor_));
}

void Site::calculateDwellTimeConstant_() {
  auto sumRates = getSumOfRates_();
  escape_time_constant_ = 1.0 / sumRates;
}

double Site::getSumOfRates_() {
  double sum = 0.0;
  for (auto rate : neighRates_) {
    sum += *(rate.second);
  }
  return sum;
}

}
