#include "common.hpp"
#include "AtomConstants.hpp"
#include "Vector3d.hpp"
#ifndef ATOM_H_
#define ATOM_H_
namespace fragdock {
  class Atom : public Vector3d {
    int atom_id;
    int xs_type;
  public:
    Atom() {}
    Atom(int atom_id, const Vector3d &pos, int xs_type) : atom_id(atom_id), xs_type(xs_type) {
      this->x = pos.x; this->y = pos.y; this->z = pos.z;
    }
    int getId() const { return atom_id; }
    void setId(int id) { atom_id = id; }
    void setPos(const Vector3d& pos) {
      x = pos.x; y = pos.y; z = pos.z;
    }
    int getXSType() const { return xs_type; }
    bool operator<(const Atom& o) const { return atom_id < o.atom_id; }
    friend std::ostream& operator<< (std::ostream& os, const Atom& atom);
  };
}
#endif
