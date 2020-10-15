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
}
}

#endif  // MYTHICAL_CONSTANTS_HPP
