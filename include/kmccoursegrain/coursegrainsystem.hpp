#ifndef KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP
#define KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP

#include <vector>
#include <memory>
#include <map>

namespace kmccoursegrain{

class KMC_Site;
class KMC_Cluster;
class Particle;

typedef std::shared_ptr<Particle> ParticlePtr;
typedef std::shared_ptr<KMC_Site> SitePtr;
typedef std::shared_ptr<KMC_Cluster> ClusterPtr;

/**
 * \brief Course Grain System allows abstraction of renormalization of sites
 *
 * This class simply uses pointers to the relevant hop rates. It will simulate a
 * particle hopping through the system of sites. If a large number of compute
 * cycles are expended moving a particle between two low energy sites they will
 * be course grained, or renormalized so that the probabilities and time spent
 * on each site will be the same but the number of compute cycles will be 
 * significantly reduced. 
 **/
class CourseGrainSystem{

  public:
  
    /**
     * \brief Constuctor for course grained system
     *
     * The constructor by default seeds the random number generator based on the
     * current time. Furthermore, course graining of two sites is by default set
     * to a threshold of 20. Meaning a particle must remember moving back and
     * forth between two sites at least 20 times before the sites are course 
     * grained. 
     **/
    CourseGrainSystem() : seed_set_(false), courseGrainingThreshold_(20) {};

    /**
     * \brief This will correctly initialize the system 
     *
     * Essentially pointers to all the rates are stored in this class. Pointers
     * are used in case any of the rates are changed, there will be no need to
     * pass them back into the couse grain object. An updated function could
     * simply be called. This function must be called before the particles can
     * be initialized `initializeParticles` and before a hopping event is called
     * on a particle `hop`. 
     *
     * \param[in] ratesOfAllSites this is a map of maps the first int is site i
     * the value of which is a second map of sites j which are all neighbors of
     * site i. The double of the final map is a pointer to the actual rate. 
     * Consider the following to understand how the maps are structured. I have
     * a site 1 which has three neighbors and a site 2 which has 2 neighbors
     *
     * site4 - site1 - site2 - site3
     *           |
     *         site5
     *
     * Each line - is composed of two rates <- and ->. So if I were to store the
     * rates in the map it would look like this
     *
     * map<int,map<int,double *>> rates;
     *
     * // Rates from site 1
     *
     * rates[1][2] = &rateFrom1to2;
     * rates[1][4] = &rateFrom1to4;
     * rates[1][5] = &rateFrom1to5;
     *
     * // Rates from site 2
     *
     * rates[2][1] = &rateFrom2to1;
     * rates[2][3] = &rateFrom2to3;
     *
     * // Rates from site 3
     *
     * rates[3][2] = &rateFrom3to2;
     *
     * // Rates from site 4
     *
     * rates[4][1] = &rateFrom4to1;
     *
     * // Rates from site 5
     *
     * rates[5][1] = &rateFrom5to1;
     *
     * Where each of the rateFrom variables is a double
     **/
    void initializeSystem(
        std::map<const int,std::map<int const,double * >> ratesOfAllSites);

    /**
     * \brief Initialize particle dwell times and future hop site id
     *
     * This function can only be called after initializing the sites. Each 
     * particle mst have also been placed on a site in the system. As in if
     * there are sites 1 2 and 3. Then the particles must exist on at least one
     * of these sites before they are passed in. The function will then update
     * their dwell times as well as providing a potential future hopping site. 
     *
     * \param[in] particles a vector of pointers to the particles
     **/
    void initializeParticles(std::vector<ParticlePtr> particles);

    /**
     * \brief Define the seed for the random number generator
     *
     * This allows the user to create reproducable results if desired. By 
     * default the seed will be determined from the time.
     *
     * \param[in] seed 
     **/
    void setRandomSeed(const unsigned long seed);

    /**
     * \brief Make the particle hop to a site in the system
     *
     * Once a particle has been initialized, i.e. it has a dwell time it is
     * located on a site in the system and it has a stored potential site it 
     * will be hopping to. This function be called on it. It will move the
     * particle if necessary and will course grain sites/renormalize sites if
     * the particle begins hopping back and forth between the sites. 
     *
     * \param[in] particle
     **/
    void hop(ParticlePtr particle);

    /**
     * \brief Threshold for course graining sites
     *
     * This function allows you to set how many times a particle will hop back
     * and forth between two sites before the sites are renormalized into a
     * cluster. By default a value of 20 is used. 
     *
     * \param[in] thershold
     **/
    void setCourseGrainThreshold(int threshold);
  private:

    /// Depicts whether a random seed has been set, to yield reproducable data
    bool seed_set_;

    /// The random seed
    unsigned long seed_;

    /// Number of times a particle moves back and forth before turned into a
    /// cluster
    int courseGrainingThreshold_;

    /// Stores smart pointers to all the sites
    std::map<int,SitePtr> sites_;    

    /// Stores smart pointers to all the clusters
    std::map<int,ClusterPtr> clusters_;

    void courseGrainSiteIfNeeded_(ParticlePtr particle);
    void mergeSiteToCluster_(const int siteId, const int clusterId);
    void createCluster_(const int siteId1,const int siteId2);
    void mergeClusters_(const int clusterId1, const int clusterId2);

};

}
#endif // KMCCOURSEGRAIN_COURSEGRAINSYSTEM_HPP
