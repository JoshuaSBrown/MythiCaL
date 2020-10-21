#ifndef MYTHICAL_CUBIC_LATTICE_HPP
#define MYTHICAL_CUBIC_LATTICE_HPP

#include "lattice.hpp"

#include <random>
#include <vector>

namespace mythical {

  namespace lattice {
   
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
          X,
          Y,
          Z
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

      private:
        int length_ = 0;
        int width_ = 0;
        int height_ = 0;
        double inter_site_distance_ = 1.0; // nm nanometers

        std::uniform_int_distribution<int> distribution_x_;
        std::uniform_int_distribution<int> distribution_y_;
        std::uniform_int_distribution<int> distribution_z_;
        std::default_random_engine generator_;

        BoundarySetting x_bound_ = BoundarySetting::Fixed;
        BoundarySetting y_bound_ = BoundarySetting::Fixed;
        BoundarySetting z_bound_ = BoundarySetting::Fixed;

        int getXPeriodic_( const int x ) const noexcept;
        int getYPeriodic_( const int y ) const noexcept;
        int getZPeriodic_( const int z ) const noexcept;

        void checkBounds_() const;
        void setDistributions_();
    };
  }
}
#endif  // MYTHICAL_CUBIC_LATTICE_HPP
