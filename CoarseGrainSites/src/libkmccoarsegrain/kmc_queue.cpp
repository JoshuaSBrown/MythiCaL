
#include "../../include/kmccoarsegrain/kmc_queue.hpp"

#include <algorithm>
#include <iostream>

using namespace std;

struct Comparitor {
    explicit Comparitor(pair<int,double> value) : value_(value) { }
    inline bool operator()(const pair<int,double> & value) const {
      return value_.second < value.second;
    }
  private:
    pair<int,double> value_;

};

namespace kmccoarsegrain {

  size_t KMC_Queue::size(){ return walker_queue_.size(); }

  void KMC_Queue::add(pair<int,double> walker){
    Comparitor comp(walker);
    auto it = find_if(walker_queue_.begin(),walker_queue_.end(),comp);
    walker_queue_.insert(it,walker);
  }

  pair<int,double> KMC_Queue::pop_current() {
    auto current = walker_queue_.front();
    walker_queue_.pop_front();
    return current;
  }

}
