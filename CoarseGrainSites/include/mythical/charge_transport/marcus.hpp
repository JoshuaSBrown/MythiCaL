#ifndef MYTHICAL_CHARGE_TRANSPORT_MARCUS_HPP
#define MYTHICAL_CHARGE_TRANSPORT_MARCUS_HPP

namespace mythical {

  namespace charge_transport {


    /**
     * @brief Marcus rate equation 
     *
     * All energies are in eV
     */
    class Marcus {
      public:

        /**
         * @brief Constructor for the Marcus Rate equation
         *
         * @param lamba - Reorganization energy
         * @param T - Temperature
         * @param H_AB - Electronic coupling
         *
         *
         * a = 2*pi/h_bar * |H_AB|^2 * 1/sqrt(4*pi*k_B*T)
         *
         * k = a * exp(-(lambda_ + DeltaG)^2/(4*pi*k_B*T));
         */
        Marcus(const double lamba, const double T, const double H_AB);

        /**
         * @brief Get the Rate
         *
         * @param E_i - The energy of the site hopping from [ eV ]
         * @param E_j - The energy of the site hopping to [ eV ]
         *
         * @return the rate k [ 1/s ]
         */
        double getRate(const double E_i, const double E_j) const noexcept;
   
      private:
        // Reorganization energy
        double lambda_;
       
        // Eponent denominator (4*pi*k_B*T) 
        double expon_denom_;

        // 2*pi/h_bar * |H_AB|^2 * 1/sqrt(4*pi*k_B*T)
        double pre_factor_;
    };
    
  }
}
#endif  // MYTHICAL_CHARGE_TRANSPORT_MARCUS_HPP
