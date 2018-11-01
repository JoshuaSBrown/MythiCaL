#include "kmc_site.hpp"
#include "../../include/kmccoursegrain/kmc_constants.hpp"
#include <chrono>
#include <set>
#include <utility>
#include <algorithm>

using namespace std;

namespace kmccoursegrain {

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
    : totalVisitFreq_(0),
      clusterId_(constants::unassignedId),
      siteOccupied_(false) {

  auto seed = chrono::system_clock::now().time_since_epoch().count();
  randomEngine_ = mt19937(seed);
  randomDistribution_ = uniform_real_distribution<double>(0.0, 1.0);
}

void KMC_Site::setRandomSeed(const unsigned long seed) {
  randomEngine_ = mt19937(seed);
}

void KMC_Site::setRatesToNeighbors(unordered_map<int, double*> neighRates) {
  for (auto neighAndRate : neighRates) {
    neighRates_[neighAndRate.first] = neighAndRate.second;
  }
  calculateDwellTimeConstant_();
  calculateProbabilityHopToNeighbors_();
}

void KMC_Site::addNeighRate(const pair<int, double*> neighRate) {

  if (neighRates_.count(neighRate.first)) {
    throw invalid_argument("That neighbor has already been added.");
  }
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
  if (neighRates_.count(neighSiteId) == 0) {
    string err = "Error the site Id " + to_string(neighSiteId) +
                 " is not a "
                 "neighbor of site " +
                 to_string(getId());
    throw invalid_argument(err);
  }
  return *(neighRates_[neighSiteId]);
}

vector<int> KMC_Site::getNeighborSiteIds() const {
  vector<int> neighborIds;
  for (auto neighId : neighRates_) neighborIds.push_back(neighId.first);
  return neighborIds;
}

double KMC_Site::getDwellTime() {
  double number = randomDistribution_(randomEngine_);
  return (-log(number) * escapeTimeConstant_);
}

int KMC_Site::pickNewSiteId() {
  double number = randomDistribution_(randomEngine_);
  double threshold = 0.0;
  for (auto pval : probabilityHopToNeighbor_) {
    threshold += pval.second;
    if (number < threshold) return pval.first;
  }
  string err =
      "Error cummulitive probability distribution is flawed or "
      " the random number generator has calculated a value greater than 1, or"
      " the site has no neighbors";

  throw invalid_argument(err);
}

double KMC_Site::getProbabilityOfHoppingToNeighboringSite(
    const int neighSiteId) 
{
  // Use neighRates because random access of map is faster than vector
  if (neighRates_.count(neighSiteId) == 0) {
    string err = "Error site " + to_string(neighSiteId) +
                 " is not a neighbor of " + to_string(getId());
    throw invalid_argument(err);
  }

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
                         const kmccoursegrain::KMC_Site& site) {
  os << "Site Id: " << site.getId() << endl;
  os << "Cluster Id: " << site.clusterId_ << endl;
  os << "Total Visit Frequency: " << site.totalVisitFreq_ << endl;
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
  escapeTimeConstant_ = 1.0 / sumRates;
}

double KMC_Site::getSumOfRates_() {
  double sum = 0.0;
  for (auto rate : neighRates_) {
    sum += *(rate.second);
  }
  return sum;
}

}
