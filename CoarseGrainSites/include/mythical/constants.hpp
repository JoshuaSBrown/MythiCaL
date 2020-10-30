#ifndef MYTHICAL_CONSTANTS_HPP
#define MYTHICAL_CONSTANTS_HPP

#include <limits>

namespace mythical {
  namespace constants {
    /// This constant is used to determine if an id has been set or not any id
    /// that has not been set will be assigned the unassignedId
    const int unassignedId = std::numeric_limits<int>::min();
    const int inf_iterations = std::numeric_limits<int>::max();
    const double unassigned_value = std::numeric_limits<double>::min();
    // Boltzmann's constant [ eV / K ]
    const double k_B = 8.617333262145e-5;
    // Square Root of Boltzmann's constant
    const double SRk_B = 0.009282959260;
    // Plank's constant [ eV * s ]
    const double h_bar = 6.582119569e-16;
    // PI
    const double PI = 3.14159265359;
    // Square Root of PI
    const double SRPI = 1.77245385091;
  }
}

#endif  // MYTHICAL_CONSTANTS_HPP
