
#include <kmccoursegrain/site.hpp>

using namespace std;
using namespace kmccoursegrain;

void Site::setRatesToNeighbors(map<int const,double&> neighRates){
  neighRates_ = neighRates;
}
