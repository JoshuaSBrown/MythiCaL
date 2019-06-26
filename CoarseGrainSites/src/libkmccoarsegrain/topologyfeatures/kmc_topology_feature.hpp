#ifndef KMCCOARSEGRAIN_KMC_TOPOLOGY_FEATURE_HPP
#define KMCCOARSEGRAIN_KMC_TOPOLOGY_FEATURE_HPP

#include <iostream>
#include <unordered_map>
#include <math.h>
#include <memory>
#include <random>
#include <vector>

#include "../../../include/kmccoarsegrain/kmc_constants.hpp"
#include "../identity.hpp"

namespace kmccoarsegrain {

/**
 * \brief TopologyFeature Class
 *
 * This class keeps track of all information related to a feature and it's
 * neighbors. It is an internal class meaning it is not meant to be used by
 * the public. 
 **/
class KMC_TopologyFeature { //: public virtual Identity {

  protected:

		int id_ = -1;
  /**
   * \brief Keeps track of the total number of time the feature has been
   * visited by a particle
   **/
  int total_visit_freq_;

  /**
   * \brief Keeps track of whether a feature is currently occupied
   *
   * A value of 0 indicates it is not occupied
   **/
  int occupied_;

  /**
   * \brief The time constant of the feature, dwell time constant
   *
   * This value is calcualted by considering the rates off the feature. If all
   * the rates off the feature are very small the time constant will be large.
   * If the rates off are fast than the time constant will be small. It has
   * an inverse relation with the rates off.
   */
  double escape_time_constant_;

  /**
   * \brief This is the random number engine
   *
   * This engine creates pseudo random numbers, it has been initialized
   * from the time.
   **/
  std::mt19937 random_engine_;

  /**
   * \brief Distribution to be used when picking random numbers
   **/
  std::uniform_real_distribution<double> random_distribution_;

  /// Create function pointer variables

  void (*occupy_ptr_)(KMC_TopologyFeature *);
  void (*occupy_siteId_ptr_)(KMC_TopologyFeature*,const int&);
 
  void (*vacate_ptr_)(KMC_TopologyFeature *);
  void (*vacate_siteId_ptr_)(KMC_TopologyFeature *,const int&);

  bool (*isOccupied_ptr_)(const KMC_TopologyFeature *);
  bool (*isOccupied_siteId_ptr_)(const KMC_TopologyFeature *,const int&);

  void (*remove_ptr_)(KMC_TopologyFeature *,const int &);

  friend void occupyTopology_(KMC_TopologyFeature*);
  friend void occupyTopology_(KMC_TopologyFeature*,const int&);

  friend void vacateTopology_(KMC_TopologyFeature*);
  friend void vacateTopology_(KMC_TopologyFeature*,const int&);

  friend bool isOccupiedTopology_(const KMC_TopologyFeature*);
  friend bool isOccupiedTopology_(const KMC_TopologyFeature*,const int&);

  friend void removeWalker_(KMC_TopologyFeature*,const int&);

 public:
  KMC_TopologyFeature();

	void setId(const int id){ id_ = id;}
	int getId() const noexcept { return id_;}
  virtual ~KMC_TopologyFeature() {};
  /**
   * \brief Set the seed for the random number generator
   *
   * By default this is initialized in the constructor from the time.
   * However, having the ability to set it allows to reproducably test the
   * class.
   *
   * \param[in] seed a random number seed
   **/
  void setRandomSeed(const unsigned long seed);

  /**
   * \brief Determine if site is occupied by a particle
   *
   * \return True if it is occupied, False if it is not occupied
   **/
  bool isOccupied() const { return isOccupied_ptr_(this); }
  bool isOccupied(const int& siteId) const 
  { return isOccupied_siteId_ptr_(this,siteId); }

  /**
   * \brief Indicate that the site is no longer occupied by a particle
   **/
  void vacate() { vacate_ptr_(this); }
  void vacate(const int& siteId) { vacate_siteId_ptr_(this,siteId); }

  /**
   * \brief Remove a random walker from the site 
   *
   * This is not the same as vacate. If a walker is removed from the system this
   * function should be called, as it ensures that any internal variables 
   * tracking the walker are appropriately reset.
   **/
  void removeWalker(const int & walker_id,const int & site_id) { 
    vacate(site_id);
    remove_ptr_(this,walker_id);
  }

  /**
   * \brief The action of occupying a site
   *
   * The site visitation frequency is incremented and the site is set to the
   * occupied status. 
   **/

  void occupy() { 
    occupy_ptr_(this); 
  }
  void occupy(const int& siteId) { 
    occupy_siteId_ptr_(this,siteId); 
  }

  /**
   * \brief Set the status of the site to occupied
   *
   * Does not increment the visit frequency. 
   **/
  void setToOccupiedStatus(){ occupied_=true; }
  void setToUnoccupiedStatus(){ occupied_=false; }

  /**
   * \brief Return the dwell time of the site
   *
   * \return A double which is the dwell time
   **/
  double getTimeConstant() const { return escape_time_constant_; }

  /**
   * \brief Return the hop time of the site
   *
   * This is different from the dwell time as it is multiplied by a random
   * number, each time this function is called it should give a different
   * time. This is where the Monte Carlo part of the simulation comes into
   * play.
   *
   * \return A time indicating how long a particle was on the site before it
   * hopped
   **/
  virtual double getDwellTime(const int & walker_id);
  //virtual double getDwellTime();

  /**
   * \brief Returns the id of a neighboring site
   *
   * Again this is based on a random number generator. A neighboring site id
   * is returned based on the rates to the neighboring site. If the rate to
   * the neighbor is large than it will be more likely that the site id
   * associated with that neighbor will be returned. Will throw an error if
   * the site has no neighbors.
   *
   * \return site id of a neigboring site
   **/
  virtual int pickNewSiteId(const int & walker_id) = 0;
  virtual int pickNewSiteId() { return -1;}

  virtual void setVisitFrequency(const int & frequency) 
  { total_visit_freq_ = frequency;}
  virtual void setVisitFrequency(const int & frequency, int) 
  { total_visit_freq_ = frequency;}

  virtual int getVisitFrequency() 
  { return total_visit_freq_; }
  virtual int getVisitFrequency(int) 
  { return total_visit_freq_; }


};

}

#endif  // KMCCOARSEGRAIN_KMC_TOPOLOGY_FEATURE_HPP
