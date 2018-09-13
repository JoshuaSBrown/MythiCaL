#include <iostream>
#include <cassert>
#include <vector>
#include <memory>

#include <kmccoursegrain/site.hpp>

using namespace std;
using namespace kmccoursegrain;

int main(void){

  cout << "Testing: Site constructor" << endl;
  {
    Site site;
  }

  cout << "Testing: Site Id setter" << endl;
  {
    Site site;
    site.setId(0);
  }
	
  cout << "Testing: Site Id getter" << endl;
  {
    Site site;
    site.setId(0);
    assert(site.getId()==0);

    bool fail = false;
    Site site2;
    try {
      site2.getId();
    }catch(...){
      fail = true;
    }
    assert(fail);
  }

	return 0;
}
