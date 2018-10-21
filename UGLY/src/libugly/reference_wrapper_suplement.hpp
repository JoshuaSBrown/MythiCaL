#ifndef UGLY_REFERENCEWRAPSUPPLEMENTAL_HPP
#define UGLY_REFERENCEWRAPSUPPLEMENTAL_HPP

#include <functional> 

namespace std {
  template<typename T>
    bool operator<(reference_wrapper<T> a, reference_wrapper<T> b){
      return a.get()<b.get();
    }

  template<typename T>                                                         
    bool operator==(
        const reference_wrapper<const T> a,
        const reference_wrapper<const T> b){
    return a.get()==b.get();
  }

  template<typename T>
    bool operator==(
        const reference_wrapper<T> a,
        const reference_wrapper<const T> b){
      return a.get()==b.get();
    }

  template<typename T>
    bool operator==(
        reference_wrapper<T> a,
        const reference_wrapper<T> b){
      return a.get()==b.get();
    }

  template<typename T>
    bool operator==( const reference_wrapper<T> a,const T b){
      return a.get()==b;
    }

}

#endif // UGLY_REFERENCEWRAPSUPPLEMENTAL_HPP
