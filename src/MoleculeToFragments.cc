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

  Molecule gen_new_mol(const Molecule& mol, const std::vector<int>& id_map) {
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
}

namespace fragdock {
  vector<Fragment> DecomposeMolecule(const Molecule &mol, 
                                     int max_ring_size, 
                                     bool merge_solitary,
                                     bool dummy_atom) {
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
    // merge solitary atoms
    // !!! depends on the order of `bonds` !!!
    for (int i = 0; i < bonds.size(); i++) {
      if (merge_solitary == false) break;

      // generate connected molecule
      const Bond &bond = bonds[i];
      int a = bond.atom_id1;
      int b = bond.atom_id2;

      // bonds connected to hydrogen atoms will not be treated
      if (atoms[a].getXSType() == XS_TYPE_H or atoms[b].getXSType() == XS_TYPE_H) continue;

      // check if already united
      if (uf.same(a, b)) continue; // no longer needed to do anything

      // try to unite atoms a and b
      // united_set: a set of atom ids in a (sub) fragment which includes a and b 
      // prev_set_a: a set of atom ids in a (sub) fragment which includes only a and not b
      // prev_set_b: a set of atom ids in a (sub) fragment which includes only b and not a
      utils::UnionFindTree uf_temp(uf);
      uf_temp.unite(a, b);
      const vector<vector<int> > prev_sets = uf.getSets();
      const vector<vector<int> > united_sets = uf_temp.getSets();
      vector<int> united_set;
      vector<int> prev_set_a, prev_set_b;
      for (int i = 0; i < united_sets.size(); i++) {
        if (exist_in(united_sets[i], a)) {
          united_set = united_sets[i];
          break;
        }
      }
      for (int i = 0; i < prev_sets.size(); i++) {
        if (exist_in(prev_sets[i], a)) {
          prev_set_a = prev_sets[i];
        }
        if (exist_in(prev_sets[i], b)) {
          prev_set_b = prev_sets[i];
        }
      }

      Molecule united_mol = gen_new_mol(mol, united_set);
      Molecule prev_mol_a = gen_new_mol(mol, prev_set_a);
      Molecule prev_mol_b = gen_new_mol(mol, prev_set_b);
      int nRings_prev = ringDetector(prev_mol_a.size(), prev_mol_a.getBonds()).size()
        + ringDetector(prev_mol_b.size(), prev_mol_b.getBonds()).size();
      int nRings_united = ringDetector(united_mol.size(), united_mol.getBonds()).size();

      // avoid new ring generation
      if (nRings_prev != nRings_united) { // there are new rings
        continue;
      }

      bool ok = true; // is it ok to merge?
      // internal rotation test
      // check the rotation invariancy by actually rotated it
      for (int j = 0; j < united_mol.getBonds().size(); j++) {
        Molecule test_united_mol;
        test_united_mol.append(united_mol);
        if (united_mol.getBonds()[j].is_rotor != true) continue;
        test_united_mol.bondRotate(j, 1); //1 rad rotation
        if (united_mol.calcRMSD(test_united_mol) >= 1e-5) {
          ok = false;
          break;
        }
      }
      if (ok) {
        uf.unite(a, b);
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

        if (dummy_atom) {
          dummys[set_id[a]].push_back(Atom(b, vec_b, XS_TYPE_DUMMY));
          dummys[set_id[b]].push_back(Atom(a, vec_a, XS_TYPE_DUMMY));
        } else {
          dummys[set_id[a]].push_back(Atom(b, vec_b, XS_TYPE_H));
          dummys[set_id[b]].push_back(Atom(a, vec_a, XS_TYPE_H));
        }
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
                                                        bool merge_solitary,
                                                        bool dummy_atom){
    std::vector<std::vector<Fragment> > ret(mols.size());
    #pragma omp parallel for // decomposition is independent process for each molecule
    for(int i=0; i<mols.size(); i++){
      ret[i] = DecomposeMolecule(mols[i], max_ring_size, merge_solitary, dummy_atom);
    }

    return ret;
  }

}
