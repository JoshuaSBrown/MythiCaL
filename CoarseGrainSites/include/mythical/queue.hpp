#ifndef MYTHICAL_QUEUE_HPP
#define MYTHICAL_QUEUE_HPP

#include "constants.hpp"

#include <cstddef>
#include <list>

namespace mythical {

/**
 **/
class Queue {
 public:
  Queue() {};

  /**
   * \brief get the walker at the front of the queue also removes from the list
   **/
  std::pair<int,double> pop_current();

  /**
   * \brief add to the queue
   **/
  void add(std::pair<int,double> walker);

  std::size_t size();
 private:
  /**
   * The interger in the pair is the id of the walker and the double is global
   * dwell time of the walker
   **/
  std::list<std::pair<int,double>> walker_queue_;
};
}
#endif  // MYTHICAL_QUEUE_HPP
