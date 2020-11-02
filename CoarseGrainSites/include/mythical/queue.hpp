#ifndef MYTHICAL_QUEUE_HPP
#define MYTHICAL_QUEUE_HPP

#include "constants.hpp"

#include <cstddef>
#include <deque>

namespace mythical {

/**
 **/
class Queue {
 public:
  Queue() : sorted_(true) {};

  /**
   * \brief get the walker at the front of the queue also removes from the list
   **/
  std::pair<int,double> pop_current();

  /**
   * \brief add to the queue, without sorting
   **/
  void add(std::pair<int,double> walker);

  /**
   * @brief The walker is added in the correct order.
   *
   * @param walker
   */
  void sortedAdd(std::pair<int,double> walker);

  bool isSorted() const noexcept;

  void sort(); 

  std::size_t size() const noexcept ;

  const std::pair<int,double> & at(int index) const;
 private:
  /**
   * The interger in the pair is the id of the walker and the double is global
   * dwell time of the walker.
   *
   * The Deque has been chosen because will primarily be popping from the front
   * which is faster than a list and not possible for a vector, also will
   * require sorting, we do not need to worry about invalidating pointers either
   * which means a list is not needed. 
   *
   **/
  std::deque<std::pair<int,double>> walker_queue_;

  bool sorted_;
};
}
#endif  // MYTHICAL_QUEUE_HPP
