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

  void Fragment::settri() { //一直線上にない3つの原子を単に選んで返すだけ
    tri[0] = tri[1] = tri[2] = -1;
    if (size() == 1) { //size(): 原子数(Molecule.hppで定義．FragmentクラスはMoleculeクラスを継承している)
      tri[0] = 0;
      return;
    }

    // koko ayashii kamo
    tri[0] = 0; //角度のついた3組が存在しない場合(O=C=Oなど)
    tri[1] = 1;
    for (int i = 0; i < size(); ++i) {
      assert(atoms[i].getId() == i);
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
    // logs::lout << tri[0] << " " << tri[1] << " " << tri[2] << std::endl;
  }

  void Fragment::settri(const Fragment& temp) {
    tri[0] = temp.tri[0];
    tri[1] = temp.tri[1];
    tri[2] = temp.tri[2];
  }

  // void Fragment::setrotid(const std::vector<Vector3d>& rots) {
  //   fltype theta, phi, psi;
  //   getNormalizeRot(psi, phi, theta);
  //   theta = -theta;
  //   phi = -phi;
  //   psi = -psi;

  //   Vector3d x(1, 0, 0);
  //   Vector3d y(0, 1, 0);
  //   Vector3d z(0, 0, 1);
  //   x.rotate(theta, phi, psi);
  //   y.rotate(theta, phi, psi);
  //   z.rotate(theta, phi, psi);
  //   fltype mi = 1e6;
  //   int minid = -1;
  //   for (int i = 0; i < rots.size(); ++i) {
  //     Vector3d rx(1, 0, 0);
  //     Vector3d ry(0, 1, 0);
  //     Vector3d rz(0, 0, 1);
  //     rx.rotate(rots[i]);
  //     ry.rotate(rots[i]);
  //     rz.rotate(rots[i]);
  //     fltype d = (x - rx).norm() + (y - ry).norm() + (z - rz).norm();
  //     if (d < mi) {
  //       mi = d;
  //       minid = i;
  //     }
  //   }
  //   // logs::lout << mi << " " << minid << std::endl;
  //   assert(minid != -1);
  //   rotid = minid;
  // }
  Vector3d Fragment::getRot() {
    fltype theta, phi, psi;
    getNormalizeRot(theta, phi, psi);
    return Vector3d(-psi, -phi, -theta);
  }
  Vector3d Fragment::getPos() const {
    assert(tri[0] != -1);
    return atoms[tri[0]];
  }

  void Fragment::getNormalizeRot(fltype& theta, fltype& phi, fltype& psi) {
    assert(tri[0] != -1);
    Vector3d mv = atoms[tri[0]];
    translate(-mv);
    assert((tri[1] == -1) == (size() == 1));
    theta = 0;
    phi = 0;
    psi = 0;
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
  }

  void Fragment::normalize() {
    fltype theta, phi, psi;
    getNormalizeRot(theta, phi, psi);
    // translate(-atoms[tri[0]]);
    translate(-getCenter());
    rotate(theta, phi, psi);
  }

}
