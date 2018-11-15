#ifndef KMCCOURSEGRAIN_KMC_CONSTANTS_HPP
#define KMCCOURSEGRAIN_KMC_CONSTANTS_HPP

#include <limits>

namespace kmccoursegrain {
namespace constants {
/// This constant is used to determine if an id has been set or not any id
/// that has not been set will be assigned the unassignedId
const int unassignedId = std::numeric_limits<int>::min();
}
}

#endif  // KMCCOURSEGRAIN_KMC_CONSTANTS_HPP
