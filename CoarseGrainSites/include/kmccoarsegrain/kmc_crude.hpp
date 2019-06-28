#ifndef KMCCOARSEGRAIN_KMC_CRUDE_HPP
#define KMCCOARSEGRAIN_KMC_CRUDE_HPP

#include <list>
#include <memory>
#include <random>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace kmccoarsegrain {

	class KMC_Walker;
	//class KMC_Dynamic_Topology;

	class KMC_Crude {

		public:
			KMC_Crude();
			~KMC_Crude();

			void initializeSystem(std::unordered_map<int, std::unordered_map<int, double>> &ratesOfAllSites);

	
			/**
			 * @brief Initialize the walker and return global times each walker will 
			 * spend on its current site. 
			 *
			 * @param walkers
			 *
			 * @return 
			 */
			std::list<std::pair<int,double>> initializeWalkers(std::vector<std::pair<int,KMC_Walker>>& walkers);

			void setRandomSeed(const unsigned long seed);


			void hop(std::pair<const int, KMC_Walker>& walker);
			void hop(KMC_Walker& walker);

			void removeWalkerFromSystem(std::pair<int,KMC_Walker>& walker);
			void removeWalkerFromSystem(KMC_Walker& walker);

			int getVisitFrequencyOfSite(int siteId);

		private:

			// Random number generator
			std::unordered_set<int> site_occupied_;
			std::mt19937 rand_num_gen_;
			std::uniform_real_distribution<double> distribution_;
			/// Must be a pointer so that we do not have to include private header
			//std::unique_ptr<KMC_Dynamic_Topology> topology_;
			
			std::unordered_map<int,int> visit_freq_;
			std::unordered_map<int,double> sojourn_times_;
			// Cummulitive probability distribution of hopping to neighbors 
			std::unordered_map<int,std::vector<std::pair<int,double>>> 
				cpd_neighbor_hop_;

			
	};
}
#endif  // KMCCOARSEGRAIN_KMC_CRUDE_HPP
