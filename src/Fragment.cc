#include "Fragment.hpp"

namespace fragdock {
  Fragment::Fragment(int id, const std::vector<Atom>& atoms) : id(id) {
    this->atoms = atoms;
    int ma = 0;
    for (auto& a : atoms) {
      if (ma < a.getId())
        ma = a.getId();
    }
    bond_ids.resize(ma + 1);
    // calcRadius();
    //this->translate(-getCenter());
  }

  void Fragment::settri() {
    tri[0] = tri[1] = tri[2] = -1;
    if (size() == 1) {
      tri[0] = 0;
      return;
    }

    // koko ayashii kamo
    tri[0] = 0;
    tri[1] = 1;
    for (int i = 0; i < size(); ++i) {
      assert(atoms[i].getId() == i); /* assumption: already renumbered */
    }

    // fltype minrad = 1e6;
    for (int i = 0; i < size(); ++i) {
      const Atom& a0 = atoms[i];
      int sz = bond_ids[i].size();
      for (int j = 0; j < sz; ++j) {
        for (int k = j + 1; k < sz; ++k) {
          Bond& b1 = bonds[bond_ids[i][j]];
          Bond& b2 = bonds[bond_ids[i][k]];
          int ai1 = b1.atom_id1 + b1.atom_id2 - i;
          int ai2 = b2.atom_id1 + b2.atom_id2 - i;
          const Atom& a1 = atoms[ai1];
          const Atom& a2 = atoms[ai2];

          Vector3d vec1 = a1 - a0;
          Vector3d vec2 = a2 - a0;
          fltype p = vec1.getAngle(vec2);
          const fltype pi = acos(-1.0);
          if (p > pi * 0.1 && p < pi * 0.9) {
            tri[0] = i;
            tri[1] = ai1;
            tri[2] = ai2;
            return;
          }
        }
      }
    }
  }

  Vector3d Fragment::getRot() {
    calculateNormalizeRot();
    return Vector3d(-psi, -phi, -theta);
  }

  void Fragment::calculateNormalizeRot() {
    if (theta != 0 || phi != 0 || psi != 0) return; // already calculated
    if (tri[0] == -1 && tri[1] == -1 && tri[2] == -1) settri(); // a triplet of atoms must be determined

    Vector3d mv = atoms[tri[0]];
    translate(-mv);
    assert((tri[1] == -1) == (size() == 1));
    if (tri[1] == -1) {}
    else if (tri[2] == -1) {
      Vector3d& vec = atoms[tri[1]];
      phi = -Vector3d(0, 1, 0).getAngle(Vector3d(0, vec.y, vec.z));
      theta = -Vector3d(1, 0, 0).getAngle(Vector3d(vec.x, sqrt(vec.y * vec.y + vec.z * vec.z), 0));
      if (vec.z < 0) phi = -phi;
      // rotate(theta, phi, 0);
    }
    else {
      Vector3d vec1 = atoms[tri[1]];
      Vector3d vec2 = vec1.cross(atoms[tri[2]]);
      psi = Vector3d(0, 1, 0).getAngle(Vector3d(vec2.x, vec2.y, 0));
      phi = Vector3d(0, 0, 1).getAngle(Vector3d(0, sqrt(vec2.x * vec2.x + vec2.y * vec2.y), vec2.z));
      if (vec2.x < 0) psi = -psi;
      vec1.rotate(0, phi, psi);
      theta = -Vector3d(1, 0, 0).getAngle(vec1);
      if (vec1.y < 0) theta = -theta;
      // rotate(theta, phi, psi);
    }
    translate(mv);
    return;
  }

  void Fragment::normalize_pose() {
    // translation
    translate(-getCenter());

    // rotation
    calculateNormalizeRot();
    rotate(theta, phi, psi);
  }

}
