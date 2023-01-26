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

    /* Determine three atoms which will be used to define normal rotation */
    void settri();
    /** 
     * Determine three atoms which will be used to define normal rotation 
     * Just copying from another Fragment object
     */
    void settri(const Fragment& temp);

    Vector3d getRot();
    void getNormalizeRot(fltype& theta, fltype& phi, fltype& psi);

    /* normalize position and rotation */
    void normalize_pose();

    void settempind(int ind) { tempind = ind; }
    int gettempind() const { return tempind; }

  };
}
#endif
