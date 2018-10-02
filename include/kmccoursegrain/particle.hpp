#ifndef KMCCOURSEGRAIN_PARTICLE_HPP
#define KMCCOURSEGRAIN_PARTICLE_HPP

#include <limits>
#include <map>
#include <list>
#include <vector>

namespace kmccoursegrain{

  namespace kmc_particle{

      /// Internal constant
      const int unassignedSiteId = std::numeric_limits<int>::min();

  }
  /**
   * \brief Particle class meant to be inherited 
   *
   * This class is designed to be inhertied by a class external to the library. It
   * provides the necessary interface required by functions internal to the
   * library. A key point of the class is to keep track of the particles movement
   * between sites. In this way the particle may be thought of having a memory. 
   * The size of the memory can also be controlled. In this way it is possible 
   * during runtime for the library to note if the particle becomes stuck moving
   * between two or more sites. The library can then course grain the low energy
   * sites and hopefully improve the performance.  
   **/
  class Particle{
    public: 
      Particle();

      /**
       * \brief Sets the size of the particles memory
       *
       * This allows you to change the memory of the particle. The larger the
       * memory the correspondingly larger computer resources. It essentially
       * allows the particle to record the number of unique sites visited. So
       * if a memory of 4 is set than a total of unique sites are recorded at a
       * time. 
       *
       * \param[in] numberSitesToRemember
       **/
      void setMemoryCapacity(unsigned int numberSitesToRemember);

      /**
       * \brief Returns the size of the memory
       **/
      unsigned int getMemoryCapacity() const;

      /** 
       * \brief Record the id of the site the particle currently resides on
       **/
      void occupySite(const int siteId);

      /** 
       * \brief Return the number of times the particle has visited the
       * currently occupied site
       *
       * \return the count
       **/
      int getVisitationFrequencyOfCurrentlyOccupiedSite();

      /**
       * \brief Get the site id and number of times the particle has visited it
       *
       * \return A vector of pairs is returned the first int in the pair is the 
       * id of the site, the second int is the count of the number of times the 
       * particle has visited the site. Element (0) is the most recently visited
       * site.
       **/
      std::vector<std::pair<int,int>> 
        getVisitationFrequenciesOfCurrentSites() const;

      /**
       * \brief Grabs the id of the site the particle currently resides on
       *
       * \return is of site
       **/
      int getIdOfSiteCurrentlyOccupying() const;

      /**
       * \brief Returns a list of the particles memory
       *
       * The memory consists of recording the number of unique sites visited
       * most recently and the corresponding number of times the sites have been
       * visited. 
       *
       * \return list containing pairs, first int is the site id, the second int
       * is a count indicating the number of times the particle has visited the
       * particular site. 
       **/
      std::list<std::pair<int,int>> getMemory();

      /**
       * \brief Record the site the particle will move to next
       *
       * \param[in] site id the particle will try to move to if it remains 
       * unoccupied.
       **/
      void setPotentialSite(const int siteId) { potentialSite_ = siteId; }

      /**
       * \brief Return the id of the site the particle will attemt to move to
       **/
      int getPotentialSite() const; 

      /**
       * \brief Get the dwell time 
       **/
      double getDwellTime() const { return dwelltime_; }

      /**
       * \brief Set the dwell time
       **/
      void setDwellTime(const double dwelltime) { dwelltime_ = dwelltime; }
    private:

      /// The next site the particle will try to move to
      int potentialSite_;

      /// The length of time the particle will remain on the current site
      double dwelltime_;

      /// List containing all the sites the particle has most recently visited
      std::list<std::pair<int,int>> memoryQueue_;

      /// Store the site id in the particles memory
      bool rememberSite_(int siteId);

      /**
       * \brief Make the site id `siteId` the most currently remembered site
       *
       * This function works by moving the site to the top of the list.
       *
       * \param[in] siteId id of the site most recently visited by the particle
       **/
      void refreshMemory_(const int siteId);

      /**
       * \brief Creates a new memory
       *
       * Because the memory is fixed in size, if a new memory is needed then the
       * oldest memory is removed from the list and a new memory is created such
       * that it appears at the top of the list. 
       * 
       * \param[in] siteId of a site that is not yet stored in the memory of the 
       * particle but has just been visited
       **/
      void createNewMemory_(const int siteId);

      std::list<std::pair<int,int>>::iterator 
        getMemoryIterator_(const int siteId);
  };

}
#endif // KMCCOURSEGRAIN_PARTICLE_HPP
