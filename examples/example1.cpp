#include <cassert> 
#include <vector>
#include <random>

using namespace std;

class Lattice{
  public:
    Lattice() {};

    Lattice(const int length, 
        const int width, 
        const int height,
        const double mean, 
        const double std_deviation
        ) {
      setNormalDistribution(mean, std_deviation);
      setDimensions(length, width, height);
      setRandomNumberGeneratorsForXPlane();
      assignEnergies();
    };

    ~Lattice(){};

    void setNormalDistribution(const double mean, const double std_deviation){
      mean_ = mean;
      std_deviation_ = std_deviation;
    }

    void setDimensions(const int length, const int width, const int height) {
      length_ = length;
      width_ = width;
      height_ = height;
      energies_.resize(length*width*height);
    }

    void setRandomNumberGeneratorsForXPlane(){
      assert(height_>0 && width_>0);
      auto distribution_y_ = std::uniform_int_distribution<int>(0,width_-1);
      auto distribution_z_ = std::uniform_int_distribution<int>(0,height_-1);
    }

    void setInterSiteDistance(const double distance){
      inter_site_distance_ = distance;
    }

    int getLength() const noexcept { return length_; }
    int getWidth() const noexcept { return width_; }
    int getHeight() const noexcept { return height_; }

    double getInterSiteDistance() const noexcept { return inter_site_distance_;}

    double getSiteEnergy(const int index) const noexcept {
      return energies_.at(index);
    }

    int getIndex(const int x, const int y, const int z) const noexcept {
      int y_p = getYPeriodic_(y);
      int z_p = getZPeriodic_(z);
      return (z_p*length_*width_)+(y_p*length_)+x;
    }

    int getIndex(const vector<int> & site_position) const noexcept {
      return getIndex(site_position.at(0),
          getYPeriodic_(site_position.at(1)),
          getZPeriodic_(site_position.at(2)));
    }

    int getRandomSiteOnXPlane(const int x_plane) {
      assert(x_plane>=0 && x_plane<length_);
      int y_lattice_position = distribution_y_(generator_);
      int z_lattice_position = distribution_z_(generator_);
      return getIndex(x_plane,y_lattice_position,z_lattice_position);
    }

    vector<int> getPosition(const int index) const noexcept {
      int local_index = index; 
      int z = index / (length_*width_);
      local_index -= (length_*width_*z);
      int y = local_index / (length_);
      int x = local_index % length_;
      return vector<int>{x,y,z};
    }

    void assignEnergies(){
      default_random_engine generator;
      normal_distribution<double> distribution(mean_,std_deviation_);
      for(size_t index = 0;index<energies_.size();++index){
        energies_.at(index) = distribution(generator);
      }
    }

  private:
    int length_ = 0;
    int width_ = 0;
    int height_ = 0;
    double mean_ = 0.0;
    double std_deviation_ = 0.0;
    double inter_site_distance_ = 1.0; // nm nanometers
    vector<double> energies_; // eV electron volts

    uniform_int_distribution<int> distribution_y_;
    uniform_int_distribution<int> distribution_z_;
    default_random_engine generator_;

    int getZPeriodic_(const int z) const noexcept {
      if(z<0){
        return z+(z/height_)*-1*height_+height_;
      }
      return z % height_;
    }

    int getYPeriodic_(const int y) const noexcept {
      if(y<0){
        return y+width_+ (y/width_)*-1*width_;
      }
      return y % width_;
    }
};

int main() {

 

  return 0;
}
