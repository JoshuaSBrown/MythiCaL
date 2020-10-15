#include <cassert>
#include <iostream>

#include "../../../include/mythical/queue.hpp"

using namespace std;
using namespace mythical;

int main(void) {

  cout << "Testing: Queue constructor" << endl;
  { Queue kmc_queue; }

  cout << "Testing: Queue size" << endl;
  {
    Queue kmc_queue;
    assert(kmc_queue.size()==0);
  }

  cout << "Testing: Queue add" << endl;
  {
    pair<int,double> pr1{ 1, 23.1};
    pair<int,double> pr2{ 3, 10.3};
    pair<int,double> pr3{ 2, 0.13};

    Queue kmc_queue;
    kmc_queue.add(pr1);
    assert(kmc_queue.size()==1);
    kmc_queue.add(pr2);
    assert(kmc_queue.size()==2);
    kmc_queue.add(pr3);
    assert(kmc_queue.size()==3);
  }

  cout << "Testing: Queue pop_current" << endl;
  {

    pair<int,double> pr1{ 1, 23.1};
    pair<int,double> pr2{ 3, 10.3};
    pair<int,double> pr3{ 2, 0.13};

    Queue kmc_queue;
    kmc_queue.add(pr1);
    kmc_queue.add(pr2);
    kmc_queue.add(pr3);

    auto pr = kmc_queue.pop_current();
    assert(pr==pr3);
    pr = kmc_queue.pop_current();
    assert(pr==pr2);
    pr = kmc_queue.pop_current();
    assert(pr==pr1);
    assert(kmc_queue.size()==0);
  }
  return 0;
}
