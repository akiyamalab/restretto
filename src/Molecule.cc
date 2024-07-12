#include <vector>
#include <queue>
#include <set>

#include "Molecule.hpp"
#include "UnionFindTree.hpp"

namespace fragdock {
  std::ostream& operator<< (std::ostream& os, const Bond& bond) {
    os << "From:" << bond.atom_id1 << " To:" << bond.atom_id2;
    return os;
  }

  Molecule::Molecule(const std::vector<Atom> &atoms, const std::string title, const std::string smiles) : atoms(atoms), heavy_num(-1), title(title), smiles(smiles) {
    identifier = title + "," + smiles;
    int ma = 0;
    for (auto& a : atoms) {
      if (ma < a.getId())
        ma = a.getId();
    }
    bond_ids.resize(ma + 1);
    // calcRadius();
  }

  void Molecule::translate(const Vector3d &vec) {
    // center.translate(vec);
    for(int i = 0; i < atoms.size(); i++)
      atoms[i].translate(vec);
  }

  void Molecule::rotate(fltype theta, fltype phi, fltype psi) {
    // center must be (0, 0, 0)
    for(int i = 0; i < atoms.size(); i++)
      atoms[i].rotate(theta, phi, psi);
  }
  void Molecule::rotate(const Vector3d& vec) {
    // center must be (0, 0, 0)
    for(int i = 0; i < atoms.size(); i++)
      atoms[i].rotate(vec);
  }

  void Molecule::axisRotate(const Vector3d &axis, fltype th) {
    for(int i = 0; i < atoms.size(); i++)
      atoms[i].axisRotate(axis, th);
  }

  void Molecule::axisRotate(const Vector3d& base, const Vector3d& axis, fltype th, const std::vector<int>& id_set) {
    Vector3d pos = getCenter();
    translate(-base);
    for(int i = 0; i < id_set.size(); i++)
      atoms[id_set[i]].axisRotate(axis, th);
    translate(pos - getCenter());
  }

  void Molecule::bondRotate(const int bond_id, fltype th) {
    if(bond_id >= bonds.size()){
      std::cerr << "[Molecule::bondRotate] invarid bond id: " << bond_id << std::endl;
      std::cerr << "bonds.size() = " << bonds.size() << std::endl;
      for (int i = 0; i < bonds.size(); ++i) {
        std::cerr << "BondId:" << i << " " << bonds[i] << std::endl;
      }
      exit(1);
    }

    // detect the bond is in a ring or not 
    utils::UnionFindTree uf((int)atoms.size());
    for(int i=0; i<bonds.size(); i++){
      if (i == bond_id) continue;
      const Bond &bond = bonds[i];
      uf.unite(bond.atom_id1, bond.atom_id2);
    }
    if(uf.getSets()[0].size() == (int)atoms.size()) return; // do nothing
    fragdock::Vector3d bond_axis = atoms[bonds[bond_id].atom_id2] - atoms[bonds[bond_id].atom_id1];
    axisRotate(atoms[bonds[bond_id].atom_id1], bond_axis, th, uf.getSets()[0]);
    //std::cout << *this << std::endl;
  }

  void Molecule::append(const Molecule &o) {
    for(int i = 0; i < o.size(); i++)
      append(o.getAtom(i));

    for(int i = 0; i < o.bonds.size(); i++)
      append(o.bonds[i]);
    
    // calcRadius();
  }

  // !!! append : 正規化された atom/fragment を append することを考えられていない !!!
  void Molecule::append(const Atom &o) {
    if (o.getId() >= bond_ids.size()) {
      bond_ids.resize(o.getId() + 1);
    }
    atoms.push_back(o);

    // calcRadius();
  }

  void Molecule::append(const Bond &bond) {
    assert(bond.atom_id1 < bond_ids.size() && bond.atom_id2 < bond_ids.size());
    int i = bonds.size();
    bonds.push_back(bond);
    bond_ids[bond.atom_id1].push_back(i);
    bond_ids[bond.atom_id2].push_back(i);
  }

  Vector3d Molecule::getCenter() const {
    Vector3d center(0.0, 0.0, 0.0);
    int c = 0;
    for (auto& a : atoms) {
      if (a.getXSType() != XS_TYPE_H && a.getXSType() != XS_TYPE_DUMMY) {
        center += a;
        ++c;
      }
    }
    return center/(fltype)c;
  }

  fltype Molecule::getRadius() const {
    Vector3d center = getCenter();
    fltype rad = 0.0;
    for (auto& a : atoms) {
      if (a.getXSType() != XS_TYPE_H && a.getXSType() != XS_TYPE_DUMMY) {
        fltype val = (a - center).abs();
        if (rad < val) rad = val;
      }
    }
    return rad;
  }

  int Molecule::heavyNum() {
    if (heavy_num < 0) {
      heavy_num = 0;
      for (int i = 0; i < atoms.size(); i++) {
        if (atoms[i].getXSType() != XS_TYPE_H)
          heavy_num++;
      }
    }
    return heavy_num;
  }

  std::vector<std::vector<int> > Molecule::getGraphDistances() const {
    int n = 0;
    for (const Atom& a : atoms) {
      n = std::max(n, a.getId() + 1);
    }
    // int n = atoms.size();
    std::vector<std::vector<int> > dist(n, std::vector<int>(n, -1));
    for (const Atom& a : atoms) {
      int frm = a.getId();
      // assert(frm < dist.size());
      dist[frm][frm] = 0;
      std::queue<int> q;
      q.push(frm);
      while (!q.empty()) {
        int p = q.front(); q.pop();
        // assert(p < bond_ids.size());
        for (auto& id : bond_ids[p]) {
          // assert(id < bonds.size());
          const Bond& b = bonds[id];
          int to = b.atom_id1 + b.atom_id2 - p;
          // assert(b.atom_id1 == p || b.atom_id2 == p);
          // assert(to < dist.size());
          if (dist[frm][to] == -1) {
            dist[frm][to] = dist[frm][p] + 1;
            q.push(to);
          }
        }
      }
    }
    return dist;
  }



  std::ostream& operator<< (std::ostream& os, const Molecule& mol){
    std::vector<Atom> atoms = mol.atoms;
    std::vector<Bond> bonds = mol.bonds;

    for(int i=0; i<atoms.size(); i++){
      Atom atom = atoms[i];
      os << atom << std::endl;
    }

    for(int i=0; i<bonds.size(); i++){
      Bond bond = bonds[i];
      os << "BondId:" << i << " " << bond << std::endl;
    }
    return os;
  }
  void Molecule::renumbering(int newsz, const std::vector<unsigned int>& vec) {
    using namespace std;
    int n = vec.size();
    assert(size() == n);
    vector<Atom> newatoms(newsz);
    vector<Bond> newbonds;
    vector<vector<int> > newbond_ids(newsz);
    int ma = 0;
    for (int i = 0; i < n; ++i) {
      if (vec[i] == 0) continue;
      assert((int)vec[i] - 1 < newsz);
      newatoms[vec[i] - 1] = atoms[i];
      ma = max(ma, atoms[i].getId());
    }
    vector<int> dict(ma + 1, -1);
    for (int i = 0; i < n; ++i) {
      if (vec[i] == 0) continue;
      dict[atoms[i].getId()] = i;
    }
    for (int i = 0; i < bonds.size(); ++i) {
      int a = bonds[i].atom_id1, b = bonds[i].atom_id2;
      if (a > ma || dict[a] == -1 || b > ma || dict[b] == -1) continue;
      assert(vec[dict[a]] > 0 && vec[dict[b]] > 0);
      newbond_ids[vec[dict[a]] - 1].push_back(newbonds.size());
      newbond_ids[vec[dict[b]] - 1].push_back(newbonds.size());
      newbonds.push_back(Bond(vec[dict[a]] - 1, vec[dict[b]] - 1, bonds[i].is_rotor));
    }
    for (int i = 0; i < newsz; ++i) {
      newatoms[i].setId(i);
    }
    atoms = newatoms;
    bonds = newbonds;
    bond_ids = newbond_ids;
  }
  bool Molecule::isRenumbered() const {
    for (int i = 0; i < size(); ++i) {
      if (atoms[i].getId() != i) return false;
    }
    return true;
  }
  void Molecule::deleteHydrogens() {
    using namespace std;
    vector<unsigned int> vec(size(), 0);
    int cnt = 0;
    for (int i = 0; i < size(); ++i) {
      if (atoms[i].getXSType() != XS_TYPE_H) {
        ++cnt;
        vec[i] = cnt;
      }
    }
    renumbering(cnt, vec);
  }
  fltype Molecule::getNrots() const {

    for (int i = 0; i < size(); ++i) {
      assert(atoms[i].getId() == i);
    }

    int n = 0;
    for (const Atom& a : atoms) {
      n = std::max(n, a.getId() + 1);
    }
    std::vector<int> adj_heavy_num(n);
    for (const Atom& a : atoms) {
      int i = a.getId();
      for (auto& id : bond_ids[i]) {
        const Bond& b = bonds[id];
        const Atom& o = atoms[b.atom_id1 + b.atom_id2 - i];
        assert(i == b.atom_id1 || i == b.atom_id2);
        if (xs_is_heavy(o.getXSType())) {
          ++adj_heavy_num[i];
        }
      }
    }

    fltype ret = 0.0;
    for (const Atom& a : atoms) {
      if (a.getXSType() == XS_TYPE_H) continue;
      int i = a.getId();
      int cnt = 0;
      for (auto& id : bond_ids[i]) {
        const Bond& b = bonds[id];
        const Atom& o = atoms[b.atom_id1 + b.atom_id2 - i];
        if (b.is_rotor && xs_is_heavy(o.getXSType()) && adj_heavy_num[o.getId()] > 1) {
          ++cnt;
        }
      }
      ret += cnt * 0.5;
    }
    return ret;
  }

  fltype Molecule::calcRMSD(const Molecule& mol) const {
    // check if the two molecules are the same
    if (getsmiles() != mol.getsmiles()) {
      std::cerr << "[Molecule::calcRMSD] different smiles" << std::endl;
      return HUGE_VAL;
    }

    // check if the two molecules have the same graph distances
    if (getGraphDistances() != mol.getGraphDistances()) {
      std::cerr << "[Molecule::calcRMSD] different graph distances" << std::endl;
      return HUGE_VAL;
    }

    // calculate RMSD
    fltype sd = 0.0;
    for (int i = 0; i < size(); ++i) {
      if (getAtom(i).getXSType() != XS_TYPE_H and getAtom(i).getXSType() != XS_TYPE_DUMMY
          and mol.getAtom(i).getXSType() != XS_TYPE_H and mol.getAtom(i).getXSType() != XS_TYPE_DUMMY) {
        sd += (getAtom(i) - mol.getAtom(i)).norm();
      }
    }
    return std::sqrt(sd/size());
  }
}
