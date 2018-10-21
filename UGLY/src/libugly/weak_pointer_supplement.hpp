#ifndef UGLY_WEAKPOINTERSUPPLEMENTAL_HPP
#define UGLY_WEAKPOINTERSUPPLEMENTAL_HPP

#include <stdexcept>
#include <memory> 

namespace std 
{

  template<typename T>
    bool operator<(weak_ptr<T> a_ptr, weak_ptr<T> b_ptr)
    {
      if( auto a=a_ptr.lock()) if (auto b=b_ptr.lock()) return *a<*b;
      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

  template<typename T>
    bool operator==( 
        const shared_ptr<const T> a_ptr, 
        const shared_ptr<const T> b_ptr)
    {
      return *a_ptr==*b_ptr;
    }

  template<typename T>
    bool operator==( 
        const weak_ptr<const T> a_ptr, 
        const shared_ptr<const T> b_ptr)
    {
      if( auto a=a_ptr.lock()) return *a==*b_ptr;
      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

  template<typename T>
    bool operator==( 
        weak_ptr<T> a_ptr, 
        shared_ptr<T> b_ptr)
    {
      if( auto a=a_ptr.lock()) return *a==*b_ptr;
      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

  template<typename T>                                                         
    bool operator==( 
        const weak_ptr<const T> a_ptr, 
        const weak_ptr<const T> b_ptr)
    {
      if( auto a=a_ptr.lock()) if(auto b=b_ptr.lock()) return *a==*b;
      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

  template<typename T>
    bool operator==(
        const weak_ptr<T> a_ptr,
        const weak_ptr<const T> b_ptr)
    {
      if( auto a=a_ptr.lock()) if(auto b=b_ptr.lock()) return *a==*b;
      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

  template<typename T>
    bool operator==(
        weak_ptr<T> a_ptr,
        const weak_ptr<T> b_ptr)
    {
      if( auto a=a_ptr.lock()) if( auto b=b_ptr.lock()) return *a==*b;

      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

  template<typename T>
    bool operator==( const weak_ptr<T> a_ptr,const T b)
    {
      if( auto a=a_ptr.lock()) return *a==b;
      throw runtime_error("Cannot convert weak pointer to shared pointer for "
          "comparison operator.");
    }

}

#endif // UGLY_WEAKPOINTERSUPPLEMENTAL_HPP
