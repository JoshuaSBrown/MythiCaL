
#include "mythical/charge_transport/cubic_lattice.hpp"

#include <string>
#include <vector>
#include <stdexcept>

using namespace std;
namespace mythical {

  namespace lattice {

    void Cubic::checkBounds_() const {
      if ( length_ < 0 ) {
        throw std::invalid_argument("Cannot initialize Cubic Lattice with "
            "a negative length");
      } 
      if ( width_ < 0 ) {
        throw std::invalid_argument("Cannot initialize Cubic Lattice with "
            "a negative width");
      }
      if ( height_ < 0 ) {
        throw std::invalid_argument("Cannot initialize Cubic Lattice with "
            "a negative height");
      }
    }

    int Cubic::getXPeriodic_( const int x ) const noexcept{
      if(x<0){
        return x+(x/length_)*-1*length_+length_;
      }
      return x % length_;
    }

    int Cubic::getYPeriodic_( const int y ) const noexcept{
      if(y<0){
        return y+(y/width_)*-1*width_+width_;
      }
      return y % width_;
    }

    int Cubic::getZPeriodic_( const int z ) const noexcept{
      if(z<0){
        return z+(z/height_)*-1*height_+height_;
      }
      return z % height_;
    }

    void Cubic::setDistributions_() {
      distribution_x_ = std::uniform_int_distribution<int>(0,length_-1);
      distribution_y_ = std::uniform_int_distribution<int>(0,width_-1);
      distribution_z_ = std::uniform_int_distribution<int>(0,height_-1);
    }

    Cubic::Cubic(const int length, const int width, const int height) : 
      length_(length), width_(width), height_(height) {
        checkBounds_(); 
        setDistributions_();
      }

    Cubic::Cubic(const int length, 
        const int width, 
        const int height, 
        const double inter_site_distance) :
      length_(length), width_(width), height_(height),
      inter_site_distance_(inter_site_distance) {
        checkBounds_(); 
        setDistributions_();
      }

    Cubic::Cubic(const int length, 
        const int width, 
        const int height, 
        const double inter_site_distance,
        const BoundarySetting x_bound,
        const BoundarySetting y_bound,
        const BoundarySetting z_bound) :
      length_(length), width_(width), height_(height),
      inter_site_distance_(inter_site_distance),
      x_bound_(x_bound), y_bound_(y_bound), z_bound_(z_bound) {
        checkBounds_();
        setDistributions_();
      }

    int Cubic::getLength() const noexcept { return length_; }
    int Cubic::getWidth() const noexcept { return width_; } 
    int Cubic::getHeight() const noexcept { return height_; }

    int Cubic::getIndex(const int x, const int y, const int z) const {
      int x_lattice_pos = x;
      int y_lattice_pos = y;
      int z_lattice_pos = z;
      if ( x_bound_ == BoundarySetting::Periodic ) {
        x_lattice_pos = getXPeriodic_(x);
      } 
      if ( y_bound_ == BoundarySetting::Periodic ) {
        y_lattice_pos = getYPeriodic_(y);
      }
      if ( z_bound_ == BoundarySetting::Periodic ) {
        z_lattice_pos = getZPeriodic_(z);
      }
      return (z_lattice_pos*length_*width_) + 
        (y_lattice_pos*length_)+x_lattice_pos;  
    }

    int Cubic::getIndex(const std::vector<int> & site_position) const {
      if ( site_position.size() != 3 ) {
        throw std::invalid_argument("Unable to get index, vector must have "
            "3 dimensions, the provided vector has " + std::to_string(site_position.size()));
      }
      return getIndex(site_position.at(0),
          site_position.at(1), site_position.at(2)); 
    }

    int Cubic::getRandomSite(const Plane plane, const int plane_index) {
      int x_lattice_pos = plane_index;
      int y_lattice_pos = plane_index;
      int z_lattice_pos = plane_index;
      if ( plane == Plane::X ) {
        if (x_lattice_pos >=0 && x_lattice_pos < length_ ) {
          std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
            "with x pos " + to_string(plane_index) + " when the x"
            "bounds are 0 to " + to_string(length_-1);
            throw std::invalid_argument(error_msg);
        }
        y_lattice_pos = distribution_y_(generator_);
        z_lattice_pos = distribution_z_(generator_);
      } else if ( plane == Plane::Y ) {
        if (y_lattice_pos >=0 && y_lattice_pos < width_ ) {
          std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
            "with y pos " + to_string(plane_index) + " when the y"
            "bounds are 0 to " + to_string(width_-1);
            throw std::invalid_argument(error_msg);
        }
        x_lattice_pos = distribution_x_(generator_);
        z_lattice_pos = distribution_z_(generator_);
      } else {
        if (z_lattice_pos >=0 && z_lattice_pos < height_ ) {
          std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
            "with z pos " + to_string(plane_index) + " when the z"
            "bounds are 0 to " + to_string(height_-1);
            throw std::invalid_argument(error_msg);
        }
        x_lattice_pos = distribution_x_(generator_);
        y_lattice_pos = distribution_z_(generator_);
      }
      return getIndex(x_lattice_pos,y_lattice_pos,z_lattice_pos);
    }

    std::vector<int> Cubic::getPosition(int index) const {
      if ( index >= length_ * width_ * height_ ) {
        std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
          "with index " + to_string(index) + " index is out of bounds"
          ", the bounds are between 0 and " + 
          to_string(length_ * width_ * height_) + "."; 
        throw std::invalid_argument(error_msg);
      }
      int z = index / (length_*width_);
      index -= (length_*width_*z);
      int y = index / (length_);
      int x = index % length_;
      return vector<int>{x,y,z};
    }

  } // lattice
} // mythical 
