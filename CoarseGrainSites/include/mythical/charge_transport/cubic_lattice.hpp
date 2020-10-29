#ifndef MYTHICAL_CHARGE_TRANSPORT_CUBIC_LATTICE_HPP
#define MYTHICAL_CHARGE_TRANSPORT_CUBIC_LATTICE_HPP

#include "boundarysettings.hpp"

#include <random>
#include <unordered_map>
#include <vector>

namespace mythical {

  namespace charge_transport {
   
    /**
     * \brief Cubic class is a support class meant to help with charge transport
     * simulations.
     *
     * This class provides a hopefully useful interface to the user for mapping
     * the positions of sites to their indices when using a cubic lattice.
     **/
    class Cubic { 

      public:

        enum class Plane {
          XY,
          YZ,
          XZ
        };

        Cubic() = default;

        Cubic(const int length, const int width, const int height);

        Cubic(const int length, 
            const int width, 
            const int height, 
            const double inter_site_distance);

        Cubic(const int length, 
            const int width, 
            const int height, 
            const double inter_site_distance,
            const BoundarySetting x_bound,
            const BoundarySetting y_bound,
            const BoundarySetting z_bound);

        ~Cubic(){};

        int getLength() const noexcept;
        int getWidth() const noexcept;
        int getHeight() const noexcept; 

        int getIndex(const int x, const int y, const int z) const;

        int getIndex(const std::vector<int> & site_position) const;

        int getRandomSite(const Plane plane, const int plane_index);

        std::vector<int> getPosition(int index) const;

        int getX(int index) const;
        int getY(int index) const;
        int getZ(int index) const;

        /**
         * @brief Get all sites that are within the cutoff distance of **index**
         *
         * @param index - the index of the site
         * @param cutoff - the distance in which a site has to be of **index** 
         * to be considered a neighbor
         *
         * @return a vector of indices that are the neighbors
         */
        std::vector<int> getNeighbors(const int index, const double cutoff) const;
        std::vector<std::pair<int,double>> getNeighborDistances(const int index, const double cutoff) const;
        std::vector<std::pair<int,double>> getNeighborDistances(const std::vector<int> lattice_pos, const double cutoff) const;

        /**
         * @brief Get Neighbor Distances within cutoff for the whole lattice
         *
         * @param cutoff
         *
         * @return 
         */
        std::unordered_map<int, std::unordered_map<int, double>> getNeighborDistances(const double cutoff) const;

        /**
         * @brief Get the distance between two sites
         *
         * @param index1
         * @param index2
         *
         * @return the distance 
         */
        double getDistance(const int index1, const int index2) const;

      private:
        int length_;
        int width_;
        int height_;
        int total_;
        double inter_site_distance_; // nm nanometers

        std::uniform_int_distribution<int> distribution_x_;
        std::uniform_int_distribution<int> distribution_y_;
        std::uniform_int_distribution<int> distribution_z_;
        std::default_random_engine generator_;

        BoundarySetting x_bound_;
        BoundarySetting y_bound_;
        BoundarySetting z_bound_;

        int getXPeriodic_( const int x ) const noexcept;
        int getYPeriodic_( const int y ) const noexcept;
        int getZPeriodic_( const int z ) const noexcept;

        void checkPosX_(const int x) const;
        void checkPosY_(const int y) const;
        void checkPosZ_(const int z) const;
        void checkBounds_() const;
        void checkIndex_(const int index) const;
        void setDistributions_();
        // Only call if the x, y and z are gauranteed to be within bounds
        int getIndex_(const int x, const int y, const int z) const noexcept;

        std::vector<int> getPos_(int index) const noexcept;
        // Only call if indices have already been checked
        double getDistance_(const int x1, const int y1, const int z1,
            const int x2, const int y2, const int z2) const noexcept;

        std::vector<std::pair<int,double>> getNeighborDistances_(const std::vector<int> lattice_pos, const double cutoff) const;


    };
  }
}
#endif  // MYTHICAL_CHARGE_TRANSPORT_CUBIC_LATTICE_HPP
