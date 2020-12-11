
#include <catch2/catch.hpp>

#include <cassert>
#include <iostream>

#include "mythical/charge_transport/cuboid_lattice.hpp"

using namespace std;
using namespace mythical;
using namespace mythical::charge_transport;

TEST_CASE("Testing: Cuboid lattice constructors and getters","[unit]") {

  GIVEN("A cubic lattice with len=1, wid=3, hei=2"){
    Cuboid lattice(1, 3, 2);
    CHECK( lattice.getLength() == 1 );
    CHECK( lattice.getWidth() == 3 );
    CHECK( lattice.getHeight() == 2 );
    CHECK( lattice.getLatticeSpacing() == Approx(1.0) );
    CHECK( not lattice.isXPeriodic() );
    CHECK( not lattice.isYPeriodic() );
    CHECK( not lattice.isZPeriodic() );
  }

  GIVEN("A cubic lattice with len=2, wid=3, hei=4, site_dist = 2.0"){
    Cuboid lattice(2, 3, 4, 2.0);
    CHECK( lattice.getLength() == 2 );
    CHECK( lattice.getWidth() == 3 );
    CHECK( lattice.getHeight() == 4 );
    CHECK( lattice.getLatticeSpacing() == Approx(2.0) );
    CHECK( not lattice.isXPeriodic() );
    CHECK( not lattice.isYPeriodic() );
    CHECK( not lattice.isZPeriodic() );
  }

  GIVEN("A cubic lattice with len=2, wid=3, hei=4, site_dist = 2.0, y and z are Periodic"){
    Cuboid lattice(2, 3, 4, 2.0, BoundarySetting::Fixed, BoundarySetting::Periodic, BoundarySetting::Periodic);
    CHECK( lattice.getLength() == 2 );
    CHECK( lattice.getWidth() == 3 );
    CHECK( lattice.getHeight() == 4 );
    CHECK( lattice.getLatticeSpacing() == Approx(2.0) );
    CHECK( not lattice.isXPeriodic() );
    CHECK( lattice.isYPeriodic() );
    CHECK( lattice.isZPeriodic() );
  }

  GIVEN("A cubic lattice of size 2, 3, 4") {
    Cuboid lattice(2,3,4);
    THEN("The first index at position 0,0,0 should be 0"){
      CHECK( lattice.getIndex(0,0,0) == 0);
    }
    THEN("The index at position 1,0,0 should be 1"){
      CHECK( lattice.getIndex(1,0,0) == 1);
    }
    THEN("The index at position 0,1,0 should be 2"){
      CHECK( lattice.getIndex(0,1,0) == 2);
    }
    THEN("The index at position 0,0,1 should be 6"){
      CHECK( lattice.getIndex(0,0,1) == 6);
    }

    CHECK_THROWS( lattice.getIndex(2,0,0) );
    CHECK_THROWS( lattice.getIndex(0,3,0) );
    CHECK_THROWS( lattice.getIndex(0,0,4) );

    std::vector<int> pos = lattice.getPosition(0);
    CHECK( pos.at(0) == 0 );
    CHECK( pos.at(1) == 0 );
    CHECK( pos.at(2) == 0 );

    CHECK_THROWS( lattice.getPosition(2*3*4));
    pos = lattice.getPosition(2*3*4-1);
    CHECK( pos.at(0) == 1 );
    CHECK( pos.at(1) == 2 );
    CHECK( pos.at(2) == 3 );

    CHECK( lattice.getX(lattice.getIndex(1,0,0)) == 1 );
    CHECK( lattice.getX(lattice.getIndex(0,0,0)) == 0 );
    CHECK( lattice.getY(lattice.getIndex(1,2,0)) == 2 );
    CHECK( lattice.getY(lattice.getIndex(0,0,1)) == 0 );
    CHECK( lattice.getZ(lattice.getIndex(1,2,0)) == 0 );
    CHECK( lattice.getZ(lattice.getIndex(0,0,3)) == 3 );
  }
}

TEST_CASE("Testing: Cuboid lattice neighbors and distances","[unit]") {
  GIVEN("A cubic lattice of size 4, 5, 7") {
    Cuboid lattice(4,5,7);
    THEN("check that the neighbors are correctly found") {
      std::vector<int> neighbors = lattice.getNeighbors(0, 1.0); 
      CHECK( neighbors.size() == 3);

      neighbors = lattice.getNeighbors(4*5*7-1, 1.0); 
      CHECK( neighbors.size() == 3);

      // Ensure that the positions of neighbors are accurate because it is
      // a corner 
      int x_at_2 = 0;
      int y_at_3 = 0;
      int z_at_5 = 0;
      for ( int neigh_ind : neighbors ) {
        auto pos = lattice.getPosition(neigh_ind);
        if( pos.at(0) == 2 ) x_at_2++;
        if( pos.at(1) == 3 ) y_at_3++;
        if( pos.at(2) == 5 ) z_at_5++;
      }
      CHECK( x_at_2 == 1 );
      CHECK( y_at_3 == 1 );
      CHECK( z_at_5 == 1 );
    }
  }

  GIVEN("A periodic cubic lattice of size 4, 5, 7") {
    Cuboid lattice(4,5,7, 1.0, BoundarySetting::Periodic, 
        BoundarySetting::Periodic, BoundarySetting::Periodic);
    THEN("check that the neighbors are correctly found") {
      std::vector<int> neighbors = lattice.getNeighbors(0, 1.0); 
      CHECK( neighbors.size() == 6);

      neighbors = lattice.getNeighbors(4*5*7-1, 1.0); 
      CHECK( neighbors.size() == 6);
    }
  }

  GIVEN("A periodic cubic lattice of size 4, 5, 7") {
    Cuboid lattice(4,5,7, 1.0, BoundarySetting::Periodic, 
        BoundarySetting::Periodic, BoundarySetting::Periodic);
    THEN("check that the neighbor distances are correctly found") {
      std::vector<std::pair<int,double>> neigh_dists = lattice.getNeighborDistances(0, 1.0); 
      CHECK( neigh_dists.size() == 6);
      
      for ( auto neigh_dist : neigh_dists ) {
        CHECK( neigh_dist.second == Approx(1.0) );
      }
    }
  }

  GIVEN("A 1D cubic lattice of size 10, 1, 1, periodic in x direction") {
    Cuboid lattice(10,1,1, 1.0, BoundarySetting::Periodic, 
        BoundarySetting::Fixed, BoundarySetting::Fixed);
    THEN("check the closest distance between index 0 and 9") {
      CHECK( lattice.getSmallestDistance(0,9) == Approx(1.0));
    }
    THEN("check the closest distance between index 0 and 4") {
      CHECK( lattice.getSmallestDistance(0,4) == Approx(4.0));
    }
  }
}
