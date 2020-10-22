#include "mythical/charge_transport/marcus.hpp"
#include "mythical/constants.hpp"

#include <cmath>

using namespace mythical::constants;

namespace mythical {
  namespace charge_transport {

    Marcus::Marcus(const double lamba, const double T, const double H_AB) :
      lambda_(lambda), 
      expon_denom_(4.0*PI*k_B*T),
      pre_factor_(2.0*PI/h_bar * std::pow(H_AB,2.0) * 1.0/(SRPI*2*SRk_B*std::pow(T,0.5))) {};

    Marcus::getRate(const double E_i, const double E_j) const noexcept {
      // DeltaG = E_j - E_i 
      return pre_factor_ * std::exp(-1.0*std::pow(lambda_ + E_j - E_i,2.0)/expon_denom_);
    }
  }
}
