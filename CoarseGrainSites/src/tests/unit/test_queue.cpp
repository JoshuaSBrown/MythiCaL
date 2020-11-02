#include <cassert>
#include <iostream>

#include "mythical/queue.hpp"

using namespace std;
using namespace mythical;

int main(void) {

  cout << "Testing: Queue constructor" << endl;
  { Queue kmc_queue; }

  cout << "Testing: Queue size" << endl;
  {
    Queue kmc_queue;
    assert(kmc_queue.size()==0);
    assert(kmc_queue.isSorted() );
  }

  cout << "Testing: Queue add" << endl;
  {
    pair<int,double> pr1{ 1, 23.1};
    pair<int,double> pr2{ 3, 10.3};
    pair<int,double> pr3{ 2, 0.13};

    Queue kmc_queue;
    kmc_queue.add(pr1);
    // Add does not change the order
    assert(kmc_queue.size()==1);
    assert(kmc_queue.at(0).first == 1 );
    assert(kmc_queue.isSorted() == false);
    kmc_queue.add(pr2);
    assert(kmc_queue.size()==2);
    assert(kmc_queue.at(1).first == 3 );
    assert(kmc_queue.isSorted() == false);
    kmc_queue.add(pr3);
    assert(kmc_queue.size()==3);
    assert(kmc_queue.at(2).first == 2 );
  }

  cout << "Testing: Queue sort" << endl;
  {
    pair<int,double> pr1{ 1, 23.1};
    pair<int,double> pr2{ 3, 10.3};
    pair<int,double> pr3{ 2, 0.13};
    // Shortest time should go at the front
    Queue kmc_queue;
    kmc_queue.add(pr1);
    kmc_queue.add(pr2);
    kmc_queue.add(pr3);
    assert(kmc_queue.isSorted() == false);
    kmc_queue.sort();
    assert(kmc_queue.isSorted());
    // Add does not change the order
    assert(kmc_queue.at(0).first == 2 );
    assert(kmc_queue.at(1).first == 3 );
    assert(kmc_queue.at(2).first == 1 );
  }

  cout << "Testing: Queue sortedAdd" << endl;
  {
    pair<int,double> pr1{ 1, 23.1};
    pair<int,double> pr2{ 3, 10.3};
    pair<int,double> pr3{ 2, 0.13};

    Queue kmc_queue;
    kmc_queue.sortedAdd(pr1);
    assert(kmc_queue.size()==1);
    assert(kmc_queue.isSorted());
    kmc_queue.sortedAdd(pr2);
    assert(kmc_queue.size()==2);
    assert(kmc_queue.isSorted());
    kmc_queue.sortedAdd(pr3);
    assert(kmc_queue.size()==3);
    assert(kmc_queue.isSorted());

    assert(kmc_queue.at(0).first == 2 );
    assert(kmc_queue.at(1).first == 3 );
    assert(kmc_queue.at(2).first == 1 );
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
    kmc_queue.sort();

    auto pr = kmc_queue.pop_current();
    assert(pr==pr3);
    assert(kmc_queue.at(0).first == 3 );
    assert(kmc_queue.at(1).first == 1 );
    pr = kmc_queue.pop_current();
    assert(pr==pr2);
    assert(kmc_queue.at(0).first == 1 );
    pr = kmc_queue.pop_current();
    assert(pr==pr1);
    assert(kmc_queue.size()==0);
  }
  return 0;
}
