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
    fltype theta = 0, phi = 0, psi = 0;
    int tempind; /* TODO: Is it needed? */
    /* Determine three atoms which will be used to define normal rotation */
    void settri();
    /* Calculate rotation status [theta, phi, psi] of this pose */
    void calculateNormalizeRot();
  public:
    Fragment(int id, const std::vector<Atom> &atoms);
    Fragment() {};
    int getId() const { return id; }
    int gettri(int i) const { return tri[i]; }

    /* Get rotation status of this pose. If it is not calculated, calculate it. */
    Vector3d getRot();

    /* normalize position and rotation */
    void normalize_pose();

    void settempind(int ind) { tempind = ind; }
    int gettempind() const { return tempind; }

  };
}
#endif
