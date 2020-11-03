
#include "mythical/charge_transport/cubic_lattice.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>

using namespace std;
namespace mythical {

  namespace charge_transport {

    /**************************************************************************
     * Internal Methods
     *************************************************************************/
    void Cubic::checkPosX_(const int x) const {
      if ( x < 0 ) {
        throw std::invalid_argument("x position of site is negative, site must"
            " be within bounds between 0 and " + to_string(length_-1));
      } 
      if ( x >= length_ ) {
        throw std::invalid_argument("x position of site is too large, site must"
            " be within bounds between 0 and " + to_string(length_-1));
      } 
    }
    void Cubic::checkPosY_(const int y) const {
      if ( y < 0 ) {
        throw std::invalid_argument("y position of site is negative, site must"
            " be within bounds between 0 and " + to_string(width_-1));
      }
      if ( y >= width_ ) {
        throw std::invalid_argument("y position of site is too large, site must"
            " be within bounds between 0 and " + to_string(width_-1));
      }
    }
    void Cubic::checkPosZ_(const int z) const {
      if ( z < 0 ) {
        throw std::invalid_argument("z position of site is negative, site must"
            " be within bounds between 0 and " + to_string(height_-1));
      }
      if ( z >= height_ ) {
        throw std::invalid_argument("z position of site is too large, site must"
            " be within bounds between 0 and " + to_string(height_-1));
      }
    }

    std::vector<int> Cubic::getPos_(int index) const noexcept {
      int z = index / (length_*width_);
      index -= (length_*width_*z);
      int y = index / (length_);
      int x = index % length_;
      return vector<int>{x,y,z};
    }

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
    
    void Cubic::checkIndex_(const int index) const {
      if ( index < 0 ) {
        throw std::invalid_argument("Cannot access index " + std::to_string(index) + 
            " out of bounds, index must be a value between 0 and " + std::to_string(total_-1));
      }
      if ( index >= total_ ) {
        throw std::invalid_argument("Cannot access index " + std::to_string(index) + 
            " out of bounds, index must be a value between 0 and " + std::to_string(total_-1));
      } 
    }

    int Cubic::getIndex_(const int x, const int y, const int z) const noexcept {
      return (z*length_*width_) + (y*length_)+x;  
    }

    double Cubic::getDistance_(const int x1, const int y1, const int z1,
        const int x2, const int y2, const int z2) const noexcept {

      int xdiff = x1-x2;
      int ydiff = y1-y2;
      int zdiff = z1-z2;

      // Square
      xdiff = xdiff*xdiff;
      ydiff = ydiff*ydiff;
      zdiff = zdiff*zdiff;
      
      return std::pow(static_cast<double>(xdiff+ydiff+zdiff), 0.5) * inter_site_distance_;
    }

    std::vector<std::pair<int,double>> Cubic::getNeighborDistances_(const std::vector<int> pos_i, const double cutoff) const {

      // Calculate the number of sites that will be within the cutoff 
      int num_sites = static_cast<int>(std::floor(cutoff/inter_site_distance_));  

      const int x = pos_i.at(0);
      const int y = pos_i.at(1);
      const int z = pos_i.at(2);
      const int index = getIndex_(x, y, z);
      const int x_min = pos_i.at(0) - num_sites;
      const int x_max = pos_i.at(0) + num_sites;
      const int y_min = pos_i.at(1) - num_sites;
      const int y_max = pos_i.at(1) + num_sites;
      const int z_min = pos_i.at(2) - num_sites;
      const int z_max = pos_i.at(2) + num_sites;

      const int max_num_neigh =2*num_sites*2*num_sites*2*num_sites;
      std::vector<std::pair<int,double>> neighbors;
      neighbors.reserve(max_num_neigh);

      for ( int x_ind = x_min; x_ind <= x_max; ++x_ind ) {
        int x_lattice_pos = x_ind;
        if ( x_bound_ == BoundarySetting::Periodic ) {
          x_lattice_pos = getXPeriodic_(x_ind);
        }       
        if ( x_lattice_pos >= 0 && x_lattice_pos < length_ ) {
          for ( int y_ind = y_min; y_ind <= y_max; ++y_ind ) {
            int y_lattice_pos = y_ind;
            if ( y_bound_ == BoundarySetting::Periodic ) {
              y_lattice_pos = getYPeriodic_(y_ind);
            }       
            if ( y_lattice_pos >= 0 && y_lattice_pos < width_ ) {
              for ( int z_ind = z_min; z_ind <= z_max; ++z_ind ){
                int z_lattice_pos = z_ind;

                if ( z_bound_ == BoundarySetting::Periodic ) {
                  z_lattice_pos = getZPeriodic_(z_ind);
                }       
                if ( z_lattice_pos >= 0 && z_lattice_pos < height_ ) {
                  double dist = getDistance_(x_ind, y_ind, z_ind, pos_i.at(0), pos_i.at(1), pos_i.at(2));
                  if ( dist <= cutoff ) {
                    int neigh_ind = getIndex_(x_lattice_pos, y_lattice_pos, z_lattice_pos);
                    if( neigh_ind != index ) {
                      neighbors.emplace_back(neigh_ind, dist);
                    }
                  }
                }
              } // for z
            }
          } // for y
        }
      } // for x
      return neighbors;
    }


    int Cubic::getXPeriodic_( const int x ) const noexcept{
      if(x<0){
        return x+length_;
      }
      return x % length_;
    }

    int Cubic::getYPeriodic_( const int y ) const noexcept{
      if(y<0){
        return y+width_;
      }
      return y % width_;
    }

    int Cubic::getZPeriodic_( const int z ) const noexcept{
      if(z<0){
        return z+height_;
      }
      return z % height_;
    }

    void Cubic::setDistributions_() {
      distribution_x_ = std::uniform_int_distribution<int>(0,length_-1);
      distribution_y_ = std::uniform_int_distribution<int>(0,width_-1);
      distribution_z_ = std::uniform_int_distribution<int>(0,height_-1);
    }

    /**************************************************************************
     * Public Methods
     *************************************************************************/
    Cubic::Cubic(const int length, const int width, const int height) : 
      length_(length), width_(width), height_(height), 
      total_(length*width*height),
      inter_site_distance_(1.0),
      x_bound_(BoundarySetting::Fixed),
      y_bound_(BoundarySetting::Fixed),
      z_bound_(BoundarySetting::Fixed)
    {
      checkBounds_(); 
      setDistributions_();
    }

    Cubic::Cubic(const int length, 
        const int width, 
        const int height, 
        const double inter_site_distance) :
      length_(length), width_(width), height_(height),
      total_(length*width*height),
      inter_site_distance_(inter_site_distance),
      x_bound_(BoundarySetting::Fixed),
      y_bound_(BoundarySetting::Fixed),
      z_bound_(BoundarySetting::Fixed)
       {
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
      total_(length*width*height),
      inter_site_distance_(inter_site_distance),
      x_bound_(x_bound), y_bound_(y_bound), z_bound_(z_bound)
       {
        checkBounds_();
        setDistributions_();
      }

    int Cubic::getLength() const noexcept { return length_; }
    int Cubic::getWidth() const noexcept { return width_; } 
    int Cubic::getHeight() const noexcept { return height_; }

    bool Cubic::isXPeriodic() const noexcept {
      return x_bound_ == BoundarySetting::Periodic; 
    }

    bool Cubic::isYPeriodic() const noexcept {
      return y_bound_ == BoundarySetting::Periodic; 
    }

    bool Cubic::isZPeriodic() const noexcept {
      return z_bound_ == BoundarySetting::Periodic; 
    }

    double Cubic::getLatticeSpacing() const noexcept { 
      return inter_site_distance_;
    }

    int Cubic::getIndex(const int x, const int y, const int z) const {
      int x_lattice_pos = x;
      int y_lattice_pos = y;
      int z_lattice_pos = z;

      if ( x_bound_ == BoundarySetting::Periodic ) {
        x_lattice_pos = getXPeriodic_(x);
      } else {
        checkPosX_(x_lattice_pos);
      }

      if ( y_bound_ == BoundarySetting::Periodic ) {
        y_lattice_pos = getYPeriodic_(y);
      } else {
        checkPosY_(y_lattice_pos);
      }

      if ( z_bound_ == BoundarySetting::Periodic ) {
        z_lattice_pos = getZPeriodic_(z);
      } else {
        checkPosZ_(z_lattice_pos);
      }

      return getIndex_(x_lattice_pos, y_lattice_pos, z_lattice_pos);  
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
      if ( plane == Plane::YZ ) {
        if (x_lattice_pos < 0 && x_lattice_pos >= length_ ) {
          std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
            "with x pos " + to_string(plane_index) + " when the x"
            "bounds are 0 to " + to_string(length_-1);
            throw std::invalid_argument(error_msg);
        }
        y_lattice_pos = distribution_y_(generator_);
        z_lattice_pos = distribution_z_(generator_);
      } else if ( plane == Plane::XZ ) {
        if (y_lattice_pos < 0 && y_lattice_pos >= width_ ) {
          std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
            "with y pos " + to_string(plane_index) + " when the y"
            "bounds are 0 to " + to_string(width_-1);
            throw std::invalid_argument(error_msg);
        }
        x_lattice_pos = distribution_x_(generator_);
        z_lattice_pos = distribution_z_(generator_);
      } else {
        if (z_lattice_pos < 0 && z_lattice_pos >= height_ ) {
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
      return getPos_(index);
    }

    int Cubic::getX(int index) const {
      if ( index >= length_ * width_ * height_ ) {
        std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
          "with index " + to_string(index) + " index is out of bounds"
          ", the bounds are between 0 and " + 
          to_string(length_ * width_ * height_) + "."; 
        throw std::invalid_argument(error_msg);
      }
      int z = index / (length_*width_);
      index -= (length_*width_*z);
      return index % length_;
    }

    int Cubic::getY(int index) const {
      if ( index >= length_ * width_ * height_ ) {
        std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
          "with index " + to_string(index) + " index is out of bounds"
          ", the bounds are between 0 and " + 
          to_string(length_ * width_ * height_) + "."; 
        throw std::invalid_argument(error_msg);
      }
      int z = index / (length_*width_);
      index -= (length_*width_*z);
      return index / (length_);
    }

    int Cubic::getZ(int index) const {
      if ( index >= length_ * width_ * height_ ) {
        std::string error_msg = "Cannot call " + std::string(__FUNCTION__) + ", "
          "with index " + to_string(index) + " index is out of bounds"
          ", the bounds are between 0 and " + 
          to_string(length_ * width_ * height_) + "."; 
        throw std::invalid_argument(error_msg);
      }
      return index / (length_*width_);
    }

    std::vector<int> Cubic::getNeighbors(const int index, const double cutoff) const {
      checkIndex_(index);
      std::vector<int> pos_i = getPosition(index);
      // Calculate the number of sites that will be within the cutoff 
      int num_sites = static_cast<int>(std::floor(cutoff/inter_site_distance_));  

      const int x = pos_i.at(0);
      const int y = pos_i.at(1);
      const int z = pos_i.at(2);
      const int x_min = x - num_sites;
      const int x_max = x + num_sites;
      const int y_min = y - num_sites;
      const int y_max = y + num_sites;
      const int z_min = z - num_sites;
      const int z_max = z + num_sites;

      const int max_num_neigh = 2*num_sites*2*num_sites*2*num_sites;
      std::vector<int> neighbors;
      neighbors.reserve(max_num_neigh);

      for ( int x_ind = x_min; x_ind <= x_max; ++x_ind ) {
        int x_lattice_pos = x_ind;
        if ( x_bound_ == BoundarySetting::Periodic ) {
          x_lattice_pos = getXPeriodic_(x_ind);
        }       
        if ( x_lattice_pos >= 0 && x_lattice_pos < length_ ) {
          for ( int y_ind = y_min; y_ind <= y_max; ++y_ind ) {
            int y_lattice_pos = y_ind;
            if ( y_bound_ == BoundarySetting::Periodic ) {
              y_lattice_pos = getYPeriodic_(y_ind);
            }       
            if ( y_lattice_pos >= 0 && y_lattice_pos < width_ ) {
              for ( int z_ind = z_min; z_ind <= z_max; ++z_ind ){
                int z_lattice_pos = z_ind;

                if ( z_bound_ == BoundarySetting::Periodic ) {
                  z_lattice_pos = getZPeriodic_(z_ind);
                }       
                if ( z_lattice_pos >= 0 && z_lattice_pos < height_ ) {
                  double dist = getDistance_(x_ind, y_ind, z_ind, pos_i.at(0), pos_i.at(1), pos_i.at(2));
                  if ( dist <= cutoff ) {
                    int neigh_ind = getIndex_(x_lattice_pos, y_lattice_pos, z_lattice_pos);
                    if( neigh_ind != index ) {
                      neighbors.push_back(neigh_ind);
                    }
                  }
                }
              } // for z
            }
          } // for y
        }
      } // for x
      return neighbors;
    }

    std::vector<std::pair<int,double>> Cubic::getNeighborDistances(const int index, const double cutoff) const {
      checkIndex_(index);
      std::vector<int> pos_i = getPosition(index);

      return getNeighborDistances_(pos_i, cutoff);
    }

    std::vector<std::pair<int,double>> Cubic::getNeighborDistances(const std::vector<int> pos, const double cutoff) const {
      if ( pos.size() != 3 ) {
        throw std::invalid_argument("Unable to get NeighborDistances, vector must have "
            "3 dimensions, the provided position vector has " + std::to_string(pos.size()));
      }
      checkPosX_(pos.at(0));
      checkPosY_(pos.at(1));
      checkPosZ_(pos.at(2));

      return getNeighborDistances_(pos, cutoff);
    }

    std::unordered_map<int, std::unordered_map<int,double>> Cubic::getNeighborDistances(const double cutoff) const {

      // Calculate the number of sites that will be within the cutoff 
      int num_sites = static_cast<int>(std::floor(cutoff/inter_site_distance_));  
      
      int total_num_sites = length_*width_*height_;
      // So as not to return redundant info, only provide distance if the 
      // neighbor has a smaller index
      std::unordered_map<int, std::unordered_map<int, double>> neigh_distances;

      for ( int index = 0; index < total_num_sites; ++index ) {
        
        std::vector<int> pos_i = getPos_(index);
        const int x = pos_i.at(0);
        const int y = pos_i.at(1);
        const int z = pos_i.at(2);

        const int x_min = x - num_sites;
        const int x_max = x;
        const int y_min = y - num_sites;
        const int y_max = y;
        const int z_min = z - num_sites;
        const int z_max = z;

        for ( int x_ind = x_min; x_ind <= x_max; ++x_ind ) {
          int x_lattice_pos = x_ind;
          if ( x_bound_ == BoundarySetting::Periodic ) {
            x_lattice_pos = getXPeriodic_(x_ind);
          }       
          if ( x_lattice_pos >= 0 && x_lattice_pos < length_ ) {
            for ( int y_ind = y_min; y_ind <= y_max; ++y_ind ) {
              int y_lattice_pos = y_ind;
              if ( y_bound_ == BoundarySetting::Periodic ) {
                y_lattice_pos = getYPeriodic_(y_ind);
              }       
              if ( y_lattice_pos >= 0 && y_lattice_pos < width_ ) {
                for ( int z_ind = z_min; z_ind <= z_max; ++z_ind ){
                  int z_lattice_pos = z_ind;

                  if ( z_bound_ == BoundarySetting::Periodic ) {
                    z_lattice_pos = getZPeriodic_(z_ind);
                  }       
                  if ( z_lattice_pos >= 0 && z_lattice_pos < height_ ) {
                    double dist = getDistance_(x_ind, y_ind, z_ind, pos_i.at(0), pos_i.at(1), pos_i.at(2));
                    if ( dist <= cutoff ) {
                      int neigh_index = getIndex_(x_lattice_pos, y_lattice_pos, z_lattice_pos);
                      if ( neigh_index != index ) {
                        neigh_distances[neigh_index][index] = dist;
                      }
                    }
                  }
                } // for z
              }
            } // for y
          }
        } // for x
      } // index
      return neigh_distances;
    }

    double Cubic::getSmallestDistance(const int index1, const int index2) const {
      checkIndex_(index1);
      checkIndex_(index2);

      std::vector<int> pos_1 = getPos_(index1);
      std::vector<int> pos_2 = getPos_(index2);
      // Move the second point if a smaller distance can be found
      const int x1 = pos_1.at(0);
      const int x2 = pos_2.at(0);
      if ( x_bound_ == BoundarySetting::Periodic ){
        const int diff1 = abs(x2 - x1 + length_);
        const int diff2 = abs(x2 - x1 - length_);
        const int diff = min(diff1, diff2); 
        pos_2.at(0) = min(diff,abs(x2 - x1));
      }
      const int y1 = pos_1.at(1);
      const int y2 = pos_2.at(1);
      if ( y_bound_ == BoundarySetting::Periodic ){
        const int diff1 = abs(y2 - y1 + width_);
        const int diff2 = abs(y2 - y1 - width_);
        const int diff = min(diff1, diff2); 
        pos_2.at(1) = min(diff,abs(y2 - y1));
      }
      const int z1 = pos_1.at(1);
      const int z2 = pos_2.at(1);
      if ( z_bound_ == BoundarySetting::Periodic ){
        const int diff1 = abs(z2 - z1 + height_);
        const int diff2 = abs(z2 - z1 - height_);
        const int diff = min(diff1, diff2); 
        pos_2.at(2) = min(diff,abs(z2 - z1));
      }
      return getDistance_(pos_1.at(0), pos_1.at(1), pos_1.at(2), pos_2.at(0), pos_2.at(1), pos_2.at(2));
    }
  } // lattice
} // mythical 
