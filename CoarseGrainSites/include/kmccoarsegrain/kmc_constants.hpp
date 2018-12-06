#ifndef KMCCOARSEGRAIN_KMC_CONSTANTS_HPP
#define KMCCOARSEGRAIN_KMC_CONSTANTS_HPP

#include <limits>

namespace kmccoarsegrain {
namespace constants {
/// This constant is used to determine if an id has been set or not any id
/// that has not been set will be assigned the unassignedId
const int unassignedId = std::numeric_limits<int>::min();
const double unassigned_value = std::numeric_limits<double>::min();
}
}

#endif  // KMCCOARSEGRAIN_KMC_CONSTANTS_HPP
