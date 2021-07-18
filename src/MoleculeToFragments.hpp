#include "common.hpp"
#include "Fragment.hpp"
#include "Molecule.hpp"
#ifndef MOLECULE_TO_FRAGMENTS_H_
#define MOLECULE_TO_FRAGMENTS_H_

namespace fragdock {
  std::vector<Fragment> DecomposeMolecule(const Molecule &mol);
  std::vector<std::vector<Fragment> > DecomposeMolecule(const std::vector<Molecule> &mols);
}
#endif
