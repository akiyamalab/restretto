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


}

namespace fragdock {
  vector<Fragment> DecomposeMolecule(const Molecule &mol, 
                                     int max_ring_size, 
                                     bool merge_solitary) { // merge_solitary is not used.
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

        // near_ids[set_id[a]].push_back(b);
        // near_ids[set_id[b]].push_back(a);

        dummys[set_id[a]].push_back(Atom(b, vec_b, XS_TYPE_DUMMY));
        dummys[set_id[b]].push_back(Atom(a, vec_a, XS_TYPE_DUMMY));
        bonds_in_frags[set_id[a]].push_back(Bond(a, b, bond.is_rotor));
        bonds_in_frags[set_id[b]].push_back(Bond(a, b, bond.is_rotor));
      }
      else {
        // bonds_in_frags[set_id[a]].push_back(bond);
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
      // for (auto& j : dummys[i]) {
      //   frag_atoms.push_back(j);
      // }
      // if (dummys[i].size() >= 2 && 2 >= (int)frag_atoms.size() - (int)dummys[i].size()) continue;

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
