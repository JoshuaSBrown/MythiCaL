#ifndef KMCCOURSEGRAIN_IDENTITY_HPP
#define KMCCOURSEGRAIN_IDENTITY_HPP

#include <stdexcept>

namespace kmccoursegrain {

  class Identity {
    private:
      int id_;
      bool id_set_;

    public:
      Identity() : id_set_(false) {}
      int getId() const {
        return (id_set_ ? id_ : throw std::runtime_error("ID not set"));
      }
      void setId( int id) { id_ = id; id_set_ = true;}
  };
}

#endif // KMCCOURSEGRAIN_IDENTITY_HPP
