
#ifndef UGLY_CONSTANTS_HPP
#define UGLY_CONSTANTS_HPP

#include <limits>

namespace ugly {
  namespace constants {
    const int unassigned = std::numeric_limits<int>::min();

    enum EdgeType {
      edge,
      undirected,
      weighted,
      directed_weighted,
      directed
    };

  }
}

#endif // UGLY_CONSTANTS_HPP
