
#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "../../libmythical/identity.hpp"

using namespace std;
using namespace mythical;

TEST_CASE("Testing: Identity","[unit]"){

  cout << "Testing: Identity constructor" << endl;
  {
    Identity ident;
  }

  cout << "Testing: Identity setter" << endl;
  {
    Identity ident;
    ident.setId(0);
  }
	
  cout << "Testing: Identity getter" << endl;
  {
    Identity ident;
    ident.setId(0);
    assert(ident.getId()==0);

    bool fail = false;
    Identity ident2;
    try {
      ident2.getId();
    }catch(...){
      fail = true;
    }
    assert(fail);
  }
}
