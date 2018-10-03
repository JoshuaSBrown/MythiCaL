
#include <cmath>
#include <random>
#include <chrono>

#include "kmc_cluster.hpp"
#include "kmc_site.hpp"
#include "log.hpp"

using namespace std;

namespace kmccoursegrain {

  /****************************************************************************
   * Constants
   ****************************************************************************/

  static int clusterIdCounter = 0;
  static int threshold = 0;

  /****************************************************************************
   * Global Functions
   ****************************************************************************/
  void setThreshold(const int n){
    if(n>0){
      threshold=n;
      LOG("Threshold is "+to_string(threshold),1);
    }else{
      throw runtime_error("ERROR in setThresh, threshold is negative");
    }
  }

  int getThreshold(){
    return threshold;
  }

  /****************************************************************************
   * Public Facing Functions
   ****************************************************************************/

  KMC_Cluster::KMC_Cluster() : Identity() {
    setId(clusterIdCounter);
    clusterIdCounter++;
    iterations_ = 3;
    resolution_ = 20;
    convergenceTolerance_ = 0.01;
    convergence_method_ = converge_by_iterations_per_site;
    auto seed = std::chrono::system_clock::now().time_since_epoch().count();
    randomEngine_ = mt19937(seed);
    randomDistribution_ = uniform_real_distribution<double>(0.0,1.0);
  }

  void KMC_Cluster::addSite(SitePtr newSite){
    if(newSite == nullptr){
      throw invalid_argument("ERROR in addSite: adding a null site");
    }
    if(sitesInCluster_.count(newSite->getId())){
      throw invalid_argument("Site has already been added to the cluster");
    }
    newSite->setClusterId(this->getId());
    sitesInCluster_[newSite->getId()]=newSite;
    
    updateProbabilitiesAndTimeConstant();
  }

  void KMC_Cluster::updateProbabilitiesAndTimeConstant(){
    solveMasterEquation_();
    calculateEscapeRatesFromSitesToTheirNeighbors_();
    calculateEscapeTimeConstant_();
  }

  vector<SitePtr> KMC_Cluster::getSitesInCluster() const {
    vector<SitePtr> sites;
    for( auto site : sitesInCluster_ ) sites.push_back(site.second);
    return sites;
  }

  double KMC_Cluster::getProbabilityOfOccupyingInternalSite(const int siteId) {
    if(!sitesInCluster_.count(siteId)){
      throw invalid_argument("the provided site is not in the cluster");
    }
    return probabilityOnSite_[siteId];
  }

  double KMC_Cluster::getTimeConstant() const {
        return escapeTimeConstant_;
  }

  void KMC_Cluster::migrateSitesFrom(ClusterPtr cluster){

    for( auto site : cluster->sitesInCluster_ ) {
      site.second->setClusterId(getId());
      sitesInCluster_[site.first] = site.second;
    }

    updateProbabilitiesAndTimeConstant();
    cluster->updateProbabilitiesAndTimeConstant();
  }

  void KMC_Cluster::setRandomSeed(const unsigned long seed){
    randomEngine_ = mt19937(seed);
  }

  int KMC_Cluster::pickNewSiteId(){
    if(hopWithinCluster_()){
      return pickInternalSite_();
    }
    return pickClusterNeighbor_();
  }

  double 
  KMC_Cluster::getProbabilityOfHoppingToNeighborOfCluster(const int neighId) {
    if(!probabilityHopToNeighbor_.count(neighId)){
      string err = "Cannot get probability of hopping to neighbor, either the"
        " site is not a neighbor or the cluster has not been converged.";
      throw invalid_argument(err);
    }
    return probabilityHopToNeighbor_[neighId];
  }

  void KMC_Cluster::setConvergenceTolerance(double tolerance){
    if(tolerance<0.0){
      throw invalid_argument("tolerance must be a positive value");
    }
    convergenceTolerance_ = tolerance;
  }

  void KMC_Cluster::setConvergenceIterations(long iterations){
    if(iterations<1){
      throw invalid_argument("number of iterations must be greater than 0.");
    }
    iterations_ = iterations;
  }

  double KMC_Cluster::getDwellTime() {
    double number = randomDistribution_(randomEngine_);
    return (-log(number)*escapeTimeConstant_/resolution_);
  }

  std::ostream& 
  operator<<(std::ostream& os, const kmccoursegrain::KMC_Cluster& cluster){
    os << "Cluster Id: "        << cluster.getId()                << endl;
    os << "Cluster visitFreq: " << cluster.visitFreqCluster_      << endl;
    os << "Number of sites in Cluster: "
                                << cluster.sitesInCluster_.size() << endl;

    os << "Sites in cluster: " << endl;
    for( auto site : cluster.sitesInCluster_ ){
      os << *(site.second) << endl;
    }
    return os;
  }

  /****************************************************************************
   * Private Internal Functions
   ****************************************************************************/

  // The firsts int is the id of site within cluster
  // second int is the id of the neighbor of site outside of cluster
  map<const int, map<const int, double>>
  KMC_Cluster::getRatesToNeighborsOfCluster_(){

    map<const int, map<const int,double>> externalRates;

    for(auto site : sitesInCluster_){
      for(auto neighId : site.second->getNeighborSiteIds()) {
        if(!siteIsInCluster(neighId)){
          externalRates[site.first][neighId]=site.second->getRateToNeighbor(neighId);
        }
      }
    } 
    return externalRates;
  }

  void KMC_Cluster::initializeProbabilityOnSites_(){
    for(auto site : sitesInCluster_ ){
      probabilityOnSite_[site.first] = \
        1.0/(static_cast<double>(sitesInCluster_.size()));
    }
    return;
  }

  void KMC_Cluster::calculateProbabilityHopToNeighbors_(){

    auto ratesToNeighbors = getRatesToNeighborsOfCluster_();

    double sumRatesOffCluster = 0.0;
    for( auto rateToNeigh : ratesToNeighbors ) {
      for( auto rate : rateToNeigh.second ){
        sumRatesOffCluster+=rate.second;
      }
    }

    auto sumDwell = 0.0;
    for( auto site : sitesInCluster_ ){
      sumDwell += site.second->getTimeConstant();
    }

    map<const int, double> probabilityHopToNeighbor;
    double total = 0.0;
    for(auto rateToNeigh : ratesToNeighbors){
      int siteHoppingFrom = rateToNeigh.first;
      for( auto rate : rateToNeigh.second ){
        int siteHoppingTo = rate.first;
        if(probabilityHopToNeighbor.count(siteHoppingTo)){
          probabilityHopToNeighbor[siteHoppingTo] +=\
            probabilityOnSite_[siteHoppingFrom]*\
            sitesInCluster_[siteHoppingFrom]->getTimeConstant()/\
            sumDwell*\
            sitesInCluster_[siteHoppingFrom]->getProbabilityOfHoppingToNeighboringSite(siteHoppingTo)/\
            sumRatesOffCluster;
        }else{

          probabilityHopToNeighbor[siteHoppingTo] =\
            probabilityOnSite_[siteHoppingFrom]*\
            sitesInCluster_[siteHoppingFrom]->getTimeConstant()/\
            sumDwell*\
            sitesInCluster_[siteHoppingFrom]->getProbabilityOfHoppingToNeighboringSite(siteHoppingTo)/\
            sumRatesOffCluster;
        }
        total += probabilityHopToNeighbor[siteHoppingTo];
      }
    }

    for(auto neighborProb : probabilityHopToNeighbor ){
      probabilityHopToNeighbor[neighborProb.first] = neighborProb.second/total;
    }
    probabilityHopToNeighbor_ = probabilityHopToNeighbor;

  }


  void KMC_Cluster::iterate_(){

    map<const int, double> temp_probabilityOnSite;

    double total = 0.0;
    for(auto site : sitesInCluster_){
      for( auto neighsite : site.second->getNeighborSiteIds()){
        if(siteIsInCluster(neighsite)){
          if(temp_probabilityOnSite.count(site.first)){
            temp_probabilityOnSite[site.first] += \
              sitesInCluster_[neighsite]->getProbabilityOfHoppingToNeighboringSite(site.first)*\
              probabilityOnSite_[neighsite];
          }else{  
            temp_probabilityOnSite[site.first] = \
              sitesInCluster_[neighsite]->getProbabilityOfHoppingToNeighboringSite(site.first)*\
              probabilityOnSite_[neighsite];
          }
          total+=temp_probabilityOnSite[site.first];
        }
      }
    }

    double total2 = 0.0;
    for(auto site : sitesInCluster_){
      probabilityOnSite_[site.first] = 
        (temp_probabilityOnSite[site.first]/total+\
         probabilityOnSite_[site.first])/2.0;
      total2 += probabilityOnSite_[site.first];
    }

    for(auto site: sitesInCluster_){
      probabilityOnSite_[site.first] = probabilityOnSite_[site.first]/total2;
    }

  }

  void KMC_Cluster::solveMasterEquation_(){

    initializeProbabilityOnSites_();

    if(convergence_method_==converge_by_iterations_per_cluster){
      for(long i = 0; i < iterations_; i++){
        iterate_(); 
      }
    }else if(convergence_method_==converge_by_iterations_per_site){

      long total_iterations = iterations_*\
        static_cast<long>(sitesInCluster_.size());

      for(long i = 0; i < total_iterations; i++){
        iterate_(); 
      }
    }else{
      double error = convergenceTolerance_*1.1;
      while(error>convergenceTolerance_){
        auto oldSiteProbs = probabilityOnSite_;
        iterate_(); 
        error = 0.0;
        for( auto site : oldSiteProbs ){
          auto diff = oldSiteProbs[site.first]-probabilityOnSite_[site.first];
          error += pow(diff,2.0);
        }
        error = pow(error,1.0/2.0);
      }
    }
    calculateProbabilityHopToNeighbors_();
  }

  map<const int,vector<pair<const int, double>>> 
  KMC_Cluster::getInternalRatesFromNeighborsComingToSite_(){

    map<const int, vector<pair<const int,double>>> internalRates;

    for(auto site : sitesInCluster_){
      for(auto neighId : site.second->getNeighborSiteIds()) {
        if(siteIsInCluster(neighId)){
          auto rateToNeigh = site.second->getRateToNeighbor(neighId); 
          pair<const int,double> rateToSite(site.first,rateToNeigh);
          internalRates[neighId].push_back(rateToSite);
        }
      }
    } 
    return internalRates;
  }

  bool KMC_Cluster::hopWithinCluster_(){
    double hop = randomDistribution_(randomEngine_);
    double inOrOutThreshold = static_cast<double>(resolution_-1)/\
                              static_cast<double>(resolution_);
    return (hop<inOrOutThreshold);
  }

  int KMC_Cluster::pickClusterNeighbor_() {

    double number = randomDistribution_(randomEngine_);
    double threshold = 0.0;
    for(auto pval : probabilityHopToNeighbor_ ){
      threshold+= pval.second;
      if(number < threshold ) return pval.first;
    }
    string err = "Cummulitive probability distribution is flawed or random "
      "number is greater than 1";
    throw invalid_argument(err);
  }

  int KMC_Cluster::pickInternalSite_(){

    double number = randomDistribution_(randomEngine_);
    double threshold = 0.0;
    for(auto pval : probabilityOnSite_ ){
      threshold+= pval.second;
      if(number < threshold ) return pval.first;
    }
    string err = "cummulitive probability distribution of sites in the cluster "
      "is flawed or random number is greater than 1";
    throw invalid_argument(err);
  }

  void KMC_Cluster::calculateEscapeRatesFromSitesToTheirNeighbors_(){

    auto ratesToNeighbors = getRatesToNeighborsOfCluster_();

    map<const int, double> temp;

    for(auto site : ratesToNeighbors ){
      for(auto neigh : site.second ){
        if(temp.count(site.first)){
          temp[site.first]=neigh.second;
        }else{
          temp[site.first]+=neigh.second;
        }
      }
    }
    escapeRateFromSiteToNeighbor_ = temp;
  }

  void KMC_Cluster::calculateEscapeTimeConstant_(){
    auto ratesToNeighbors = getRatesToNeighborsOfCluster_(); 
    escapeTimeConstant_ = 0.0;
    for(auto site : ratesToNeighbors ){
      for(auto neigh : site.second ){
        escapeTimeConstant_ += neigh.second*escapeRateFromSiteToNeighbor_[site.first];
      }
    }
  }

}
