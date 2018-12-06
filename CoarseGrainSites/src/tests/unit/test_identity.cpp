#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include "../../libkmccoarsegrain/identity.hpp"

using namespace std;
using namespace kmccoarsegrain;

int main(void){

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

	return 0;
}
