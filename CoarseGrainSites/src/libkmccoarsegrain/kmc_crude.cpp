
#include <algorithm>
#include <cassert>
#include <iostream>
#include "../../include/kmccoarsegrain/kmc_crude.hpp"
#include "../../include/kmccoarsegrain/kmc_walker.hpp"
#include "../../include/kmccoarsegrain/kmc_constants.hpp"

using namespace std;

namespace kmccoarsegrain {

	bool compareSecondItemOfPair(const pair<int,double> &x, const pair<int,double> & y){
		return x.second<y.second;                                                     
	}  

	bool sortbysec(const pair<int,double> &a,                                       
			const pair<int,double> &b)                                        
	{                                                                               
		return (a.second < b.second);                                               
	} 

	KMC_Crude::KMC_Crude() : distribution_(0.0,1.0) {}
	KMC_Crude::~KMC_Crude() {}

  const std::unordered_map<int,std::unordered_map<int,double>>& KMC_Crude::rates(){
    return *rates_;
  }

	void KMC_Crude::initializeSystem(std::unordered_map<int, std::unordered_map<int, double>> &ratesOfAllSites) {

		rates_ = &(ratesOfAllSites);
		// Find the sum of all the rates for each site, important for normalizing
		// And Calculate the sojourn times
		unordered_map<int,double> sum_rates;
		for(const pair<int,unordered_map<int,double>> & site : ratesOfAllSites){
			sum_rates[site.first] = 0.0;
			for(const pair<int,double> & neigh_and_rates : site.second){
				sum_rates[site.first]+=neigh_and_rates.second;	
			}
			time_constants_[site.first]=1.0/sum_rates.at(site.first);
		}

		for(const pair<int,unordered_map<int,double>> & site : ratesOfAllSites){
			vector<pair<int,double>> probability;
			for(const pair<int,double> & neigh_and_rates : site.second){
				probability.push_back(
						pair<int,double>(neigh_and_rates.first,
							neigh_and_rates.second/sum_rates.at(site.first)));
			}	
			sort(probability.begin(),probability.end(),sortbysec);
			// Cummulitive Probability Distribution CPD
			vector<pair<int,double>> cpd;
			double pval = 0.0;                                                  
			double value = 0.0;                                                 
			for( pair<int,double> prob : probability ){                         
				prob.second+=pval;                                                
				pval = prob.second;                                               
				cpd.push_back(pair<int,double>(prob.first,value));
				value = prob.second;                                              
			}                                                                   
			cpd_neighbor_hop_[site.first] = cpd;
			assert(cpd_neighbor_hop_[site.first].size()!=0);  	
		}

	}

	list<pair<int,double>> KMC_Crude::initializeWalkers(std::vector<std::pair<int,KMC_Walker>>& walkers){

		list<pair<int,double>> walker_global_times;
		for( pair<int,KMC_Walker> & id_and_walker : walkers ){
			int siteId = id_and_walker.second.getIdOfSiteCurrentlyOccupying();
			assert(siteId!=constants::unassignedId && "Cannot initialize walkers until they have been assigned a valid ID");
			walker_global_times.push_back(pair<int,double>(id_and_walker.first,time_constants_[siteId]*log(distribution_(rand_num_gen_)*-1.0)));	
			assert(site_occupied_.count(siteId)==0 && "Cannot insert two walkers onto the same site");
			site_occupied_.insert(siteId);
			visit_freq_[siteId] = 1;

			double random_number = distribution_(rand_num_gen_);

			for(auto it = cpd_neighbor_hop_[siteId].rbegin();
					it!=cpd_neighbor_hop_[siteId].rend();++it){

				if(random_number > it->second){
					int neighId = it->first;	
					id_and_walker.second.setPotentialSite(neighId);
					id_and_walker.second.setDwellTime(time_constants_[siteId]*log(distribution_(rand_num_gen_))*-1.0);
					break;
				}
			}

		}
		walker_global_times.sort(compareSecondItemOfPair);
		return walker_global_times;
	}

  void KMC_Crude::vacateSite(const int siteId){
    assert(site_occupied_.count(siteId) == 0 && "Cannot vacate site is empty");
    site_occupied_.erase(siteId);
  }

  void KMC_Crude::occupySite(const int siteId){
    assert(site_occupied_.count(siteId)==0 && "Cannot occupy site already occupied");
    site_occupied_.insert(siteId);
  }
	void KMC_Crude::setRandomSeed(const unsigned long seed){
		rand_num_gen_ = mt19937(seed);
	}

	void KMC_Crude::hop(std::pair<const int, KMC_Walker>& walker){
		hop(walker.second);
	}

	void KMC_Crude::hop(KMC_Walker& walker){
		double random_number = distribution_(rand_num_gen_);
		int pot_siteId = walker.getPotentialSite();
		int current_siteId =walker.getIdOfSiteCurrentlyOccupying();
		if(site_occupied_.count(pot_siteId)){
			pot_siteId = current_siteId;
		}else{
			site_occupied_.erase(current_siteId);
		}
		for(auto it = cpd_neighbor_hop_[pot_siteId].rbegin();
				it!=cpd_neighbor_hop_[pot_siteId].rend();++it){

			if(random_number > it->second){
				int neighId = it->first;	
				walker.setPotentialSite(neighId);
				walker.occupySite(pot_siteId);
				walker.setDwellTime(time_constants_[pot_siteId]*log(distribution_(rand_num_gen_))*-1.0);
				site_occupied_.insert(pot_siteId);
				if(visit_freq_.count(pot_siteId)){
					++visit_freq_[pot_siteId];
				}else{
					visit_freq_[pot_siteId] =1;
				}
				break;
			}
		}
	}
/*
	void KMC_Crude::hop_occupied(KMC_Walker& walker){
		double random_number = distribution_(rand_num_gen_);
		int current_siteId =walker.getIdOfSiteCurrentlyOccupying();

		pot_siteId = current_siteId;
		for(auto it = cpd_neighbor_hop_[pot_siteId].rbegin();
				it!=cpd_neighbor_hop_[pot_siteId].rend();++it){

			if(random_number > it->second){
				int neighId = it->first;	
				walker.setPotentialSite(neighId);
				walker.occupySite(pot_siteId);
				walker.setDwellTime(time_constants_[pot_siteId]*log(distribution_(rand_num_gen_))*-1.0);
				site_occupied_.insert(pot_siteId);
				if(visit_freq_.count(pot_siteId)){
					++visit_freq_[pot_siteId];
				}else{
					visit_freq_[pot_siteId] =1;
				}
				break;
			}
		}
	}*/
	void KMC_Crude::removeWalkerFromSystem(std::pair<int,KMC_Walker>& walker){
		removeWalkerFromSystem(walker.second);
	}
	void KMC_Crude::removeWalkerFromSystem(KMC_Walker& walker){
		site_occupied_.erase(walker.getIdOfSiteCurrentlyOccupying());
	}

	int KMC_Crude::getVisitFrequencyOfSite(int siteId){
		if(visit_freq_.count(siteId)){
			return visit_freq_.at(siteId);
		}
		return 0; 
	}

}
