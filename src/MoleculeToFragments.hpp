#include "common.hpp"
#include "Fragment.hpp"
#include "Molecule.hpp"
#include "UnionFindTree.hpp"
#ifndef MOLECULE_TO_FRAGMENTS_H_
#define MOLECULE_TO_FRAGMENTS_H_

namespace fragdock {
  bool IsMergeable(const Molecule &mol, utils::UnionFindTree uf, int a, int b);
  std::vector<Fragment> DecomposeMolecule(const Molecule &mol, 
                                          int max_ring_size=-1, 
                                          bool merge_solitary=false,
                                          bool dummy_atom=true);
  std::vector<std::vector<Fragment> > DecomposeMolecule(const std::vector<Molecule> &mols, 
                                                        int max_ring_size=-1, 
                                                        bool merge_solitary=false,
                                                        bool dummy_atom=true);
}
#endif
