#include <sstream>

#include "MoleculeToFragments.hpp"
#include "UnionFindTree.hpp"
namespace {
  using namespace std;
  using namespace fragdock;

  void dfs(int now, int start,
     const vector<vector<int> > &edges,
     vector<int> &ring,
     vector<int> &route,
     vector<bool> &done) {
    if(done[now]) {
      if(now == start and route.size() > 2) {
  ring = route;
      }
      return;
    }
    done[now] = true;

    route.push_back(now);
    for(int i = 0; i < edges[now].size(); i++)
      dfs(edges[now][i], start, edges, ring, route, done);

    route.pop_back();
  }

  vector<vector<int> > ringDetector(int atom_num, const vector<Bond> &bonds) {
    vector<vector<int> > edges = vector<vector<int> >(atom_num, vector<int>());
    for(int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      edges[a].push_back(b);
      edges[b].push_back(a);
    }

    vector<vector<int> > ret;
    for(int i = 0; i < atom_num; i++) {
      vector<bool> done = vector<bool>(atom_num, false);
      vector<int> ring = vector<int>();
      vector<int> route = vector<int>();
      dfs(i, i, edges, ring, route, done);
      if(ring.size() > 0) ret.push_back(ring);
    }

    return ret;
  }

  template<typename Container>
  bool exist_in(const Container& c, const typename Container::value_type& v) {
    return ( c.end() != std::find(c.begin(), c.end(), v) );
  }

  Molecule extract_substructure(const Molecule& mol, const std::vector<int>& id_map) {
    const vector<Atom> &atoms = mol.getAtoms();
    const vector<Bond> &bonds = mol.getBonds();

    Molecule temp_mol;
    for (int j = 0; j < id_map.size(); j++) {
      temp_mol.append(atoms[id_map[j]]);
    }
    for (int j = 0; j < bonds.size(); j++) {
      if (exist_in(id_map, bonds[j].atom_id1) && exist_in(id_map, bonds[j].atom_id2)) {
        Bond new_bond = bonds[j];
        new_bond.atom_id1 = std::find(id_map.begin(), id_map.end(), bonds[j].atom_id1) - id_map.begin();
        new_bond.atom_id2 = std::find(id_map.begin(), id_map.end(), bonds[j].atom_id2) - id_map.begin();
        temp_mol.append(new_bond);
      }
    }
    return temp_mol;
  }

  bool has_new_ring(const Molecule &mol, const vector<int> &atomids_subst_a, const vector<int> &atomids_subst_b) {
    vector<int> atomids_subst_united;
    atomids_subst_united.insert(atomids_subst_united.end(), atomids_subst_a.begin(), atomids_subst_a.end());
    atomids_subst_united.insert(atomids_subst_united.end(), atomids_subst_b.begin(), atomids_subst_b.end());

    // generate molecules that are united and previous
    Molecule united_mol = extract_substructure(mol, atomids_subst_united);
    Molecule prev_mol_a = extract_substructure(mol, atomids_subst_a);
    Molecule prev_mol_b = extract_substructure(mol, atomids_subst_b);

    int nRings_prev = ringDetector(prev_mol_a.size(), prev_mol_a.getBonds()).size()
      + ringDetector(prev_mol_b.size(), prev_mol_b.getBonds()).size();
    int nRings_united = ringDetector(united_mol.size(), united_mol.getBonds()).size();
    return nRings_prev != nRings_united;
  }

  /**
   * Calculate the angle between three atoms.
   * @param[in] a Atom a (end)
   * @param[in] b Atom b (vertex)
   * @param[in] c Atom c (end)
   * @return angle between b-a and b-c in radian [0, pi]
   */
  fltype calc_angle(const Atom &a, const Atom &b, const Atom &c) {
    Vector3d vec_ba = a - b;
    Vector3d vec_bc = c - b;

    return vec_ba.getAngle(vec_bc);
  }

  fltype calc_max_angle(const Molecule &mol, int a, int b, const vector<int> & atomids_subst_b) {
    const vector<Atom> &atoms = mol.getAtoms();
    fltype max_angle = 0.0;

    for (int j = 0; j < atomids_subst_b.size(); j++) {
      const Atom &atom_bj = atoms[atomids_subst_b[j]];
      if (atom_bj.getXSType() == XS_TYPE_H) continue;
      fltype angle = calc_angle(atoms[b], atoms[a], atom_bj);
      max_angle = max(max_angle, angle);
    }

    return max_angle;
  }

  bool is_mergeable(const Molecule &mol, const vector<int> &atomids_subst_a, const vector<int> &atomids_subst_b) {
    // avoid new ring generation
    if (has_new_ring(mol, atomids_subst_a, atomids_subst_b)) {
      return false;
    }

    // check the angle invariancy
    for (int j = 0; j < mol.getBonds().size(); j++) {
      if (!mol.getBonds()[j].is_rotor) {
        continue;
      }

      int id1 = mol.getBonds()[j].atom_id1;
      int id2 = mol.getBonds()[j].atom_id2;

      fltype max_angle_1, max_angle_2;
      if (exist_in(atomids_subst_a, id1) && exist_in(atomids_subst_b, id2)) {
        // id1 in a, id2 in b
        max_angle_1 = calc_max_angle(mol, id1, id2, atomids_subst_b); // id1-id2 vs id1-b
        max_angle_2 = calc_max_angle(mol, id2, id1, atomids_subst_a); // id2-id1 vs id2-a
      } else if (exist_in(atomids_subst_a, id2) && exist_in(atomids_subst_b, id1)) {
        // id2 in a, id1 in b
        max_angle_1 = calc_max_angle(mol, id2, id1, atomids_subst_b); // id2-id1 vs id2-b
        max_angle_2 = calc_max_angle(mol, id1, id2, atomids_subst_a); // id1-id2 vs id1-a
      } else {
        continue;
      }
      if (max_angle_1 > RAD_EPS && max_angle_2 > RAD_EPS) {
        return false;
      }
    }

    return true;
  }
}

namespace fragdock {
  vector<Fragment> DecomposeMolecule(const Molecule &mol, 
                                     int max_ring_size, 
                                     bool merge_solitary) {
    const vector<Atom> &atoms = mol.getAtoms();
    const vector<Bond> &bonds = mol.getBonds();
    // unite two atoms if the bond between them is NOT 'Single'
    utils::UnionFindTree uf((int)atoms.size());
    for (int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      if (!bond.is_rotor) uf.unite(a, b);
    }
    // unite members of ring systems
    vector<vector<int> > rings = ringDetector(atoms.size(), bonds);
    for (int i = 0; i < rings.size(); i++) {
      if (max_ring_size != -1 and rings[i].size() > max_ring_size) continue;
      for (int j = 0; j < rings[i].size(); j++)
        uf.unite(rings[i][0], rings[i][j]);
    }

    // count adjacent H
    vector<int> h_cnt(atoms.size(), false);
    for (int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      if (atoms[a].getXSType() != XS_TYPE_H)
        h_cnt[b]++;
      if (atoms[b].getXSType() != XS_TYPE_H)
        h_cnt[a]++;
    }
    // unite one atom except H
    vector<bool> done(atoms.size(), false);
    for (int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      if (atoms[a].getXSType() != XS_TYPE_H and atoms[b].getXSType() != XS_TYPE_H) {
        for (int _ = 0; _ < 2; _++) {
          if ((uf.getSize(a) == 1 and h_cnt[a] <= 2 and uf.getSize(b) > 1 and !done[b]) or
              (uf.getSize(a) == 1 and uf.getSize(b) == 1 and h_cnt[a] <= 2 and h_cnt[b] <= 2)) {
            uf.unite(a, b);
            done[a] = true;
          }
          swap(a, b);
        }
      }
    }
    // merge solitary atoms except H
    // !!! depends on the order of `bonds` !!!
    for (int i = 0; i < bonds.size(); i++) {
      if (merge_solitary == false) break;

      // generate connected molecule
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      
      if (done[a] || done[b]) continue;
      if (atoms[a].getXSType() == XS_TYPE_H || atoms[b].getXSType() == XS_TYPE_H) continue;
      if (uf.same(a, b)) continue;

      const vector<vector<int> > uf_sets = uf.getSets();
      vector<int> atomids_subst_a, atomids_subst_b;
      for (int i = 0; i < uf_sets.size(); i++) {
        if (exist_in(uf_sets[i], a)) {
          // atomID set of fragments that contain atomID a
          atomids_subst_a = uf_sets[i];
        }
        if (exist_in(uf_sets[i], b)) {
          // atomID set of fragments that contain atomID b
          atomids_subst_b = uf_sets[i];
        }
      }

      if (is_mergeable(mol, atomids_subst_a, atomids_subst_b)) {
        uf.unite(a, b);
        done[a] = true;
        done[b] = true;
      }
    }

    // unite H
    for (int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      if (atoms[a].getXSType() == XS_TYPE_H or atoms[b].getXSType() == XS_TYPE_H) uf.unite(a, b);
    }
    // get fragments
    vector<vector<int> > id_sets = uf.getSets();
    id_sets = uf.getSets();
    // atom_id to set_id
    vector<int> set_id(atoms.size());
    for (int i = 0; i < id_sets.size(); i++)
      for (int j = 0; j < id_sets[i].size(); j++)
        set_id[id_sets[i][j]] = i;

    // signle bond count (to find rotatable edge)
    vector<vector<Vector3d> > axises(atoms.size());
    for (int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      if (bond.is_rotor and
          atoms[a].getXSType() != XS_TYPE_H and
          atoms[b].getXSType() != XS_TYPE_H and
          set_id[a] == set_id[b]) {
        Vector3d vec_a = atoms[b] - atoms[a];
        Vector3d vec_b = atoms[a] - atoms[b];
        vec_a /= vec_a.abs();
        vec_b /= vec_b.abs();
        axises[a].push_back(vec_a);
        axises[b].push_back(vec_b);
      }
    }
    // create edges
    // and add dummy
    vector<vector<Atom> > dummys(id_sets.size());
    vector<vector<Bond> > bonds_in_frags(id_sets.size());

    // vector<vector<int> > near_ids(id_sets.size());
    for(int i = 0; i < bonds.size(); i++) {
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;
      if(atoms[a].getXSType() == XS_TYPE_H || atoms[b].getXSType() == XS_TYPE_H) continue;

      if(set_id[a] != set_id[b]) {
        Vector3d vec_a = atoms[a];
        Vector3d vec_b = atoms[b];

        dummys[set_id[a]].push_back(Atom(b, vec_b, XS_TYPE_DUMMY));
        dummys[set_id[b]].push_back(Atom(a, vec_a, XS_TYPE_DUMMY));
        bonds_in_frags[set_id[a]].push_back(Bond(a, b, bond.is_rotor));
        bonds_in_frags[set_id[b]].push_back(Bond(a, b, bond.is_rotor));
      }
      else {
        bonds_in_frags[set_id[a]].push_back(Bond(a, b, bond.is_rotor));
      }
    }
    // create fragments
    vector<Fragment> fragments;
    for(int i = 0; i < id_sets.size(); i++) {
      vector<Atom> frag_atoms;
      sort(dummys[i].begin(), dummys[i].end());
      int di = 0;
      for(int j = 0; j < id_sets[i].size(); j++) {
        if (atoms[id_sets[i][j]].getXSType() == XS_TYPE_H) continue;
        while (di < dummys[i].size() && dummys[i][di] < atoms[id_sets[i][j]]) {
          frag_atoms.push_back(dummys[i][di]);
          ++di;
        }
        frag_atoms.push_back(atoms[id_sets[i][j]]);
      }
      while (di < dummys[i].size()) {
        frag_atoms.push_back(dummys[i][di]);
        ++di;
      }

      Fragment frag(i, frag_atoms);
      for (auto b : bonds_in_frags[i]) frag.append(b);

      fragments.push_back(frag);
    }

    return fragments;
  }

  std::vector<std::vector<Fragment> > DecomposeMolecule(const std::vector<Molecule> &mols, 
                                                        int max_ring_size, 
                                                        bool merge_solitary){
    std::vector<std::vector<Fragment> > ret(mols.size());
    #pragma omp parallel for // decomposition is independent process for each molecule
    for(int i=0; i<mols.size(); i++){
      ret[i] = DecomposeMolecule(mols[i], max_ring_size, merge_solitary);
    }

    return ret;
  }

}
