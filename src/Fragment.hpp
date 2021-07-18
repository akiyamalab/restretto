#include "common.hpp"
#include "Molecule.hpp"
#include "Atom.hpp"
#include <iostream>
#ifndef FRAGMENT_H_
#define FRAGMENT_H_
namespace fragdock {
  class Fragment : public Molecule {
    int id;
    int tri[3] = { -1, -1, -1 };
    // int rotid;
    int tempind;
  public:
    Fragment(int id, const std::vector<Atom> &atoms);
    Fragment() {};
    int getId() const { return id; }
    int gettri(int i) const { return tri[i]; }
    void settri();
    void settri(const Fragment& temp);
    // void setrotid(const std::vector<Vector3d>& rots);
    // int getrotid() const { return rotid; }
    // void getRot(fltype& theta, fltype& phi, fltype& psi);
    Vector3d getRot();
    Vector3d getPos() const;
    // Vector3d getCenter() const { return getPos(); }
    // Vector3d getNormalizeRot();
    void getNormalizeRot(fltype& theta, fltype& phi, fltype& psi);
    void normalize();

    void settempind(int ind) { tempind = ind; }
    int gettempind() const { return tempind; }

  };
}
#endif
