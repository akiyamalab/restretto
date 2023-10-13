#include "common.hpp"
#include "AtomInterEnergyGrid.hpp"
#include "Molecule.hpp"
#include "utils.hpp"
#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_
namespace fragdock {
  class Optimizer {
    const Molecule& receptor;
  public:
    explicit Optimizer(const Molecule& receptor) : receptor(receptor) {}
    fltype optimize(Molecule &mol, const EnergyCalculator& ec) const;
  };

  class Optimizer_Grid {
    const std::vector<fragdock::AtomInterEnergyGrid>& atom_grids;
    // const Molecule& receptor;
    fltype calcInterEnergy(const Molecule &mol) const {
      fltype ret = 0.0;
      for (auto& a : mol.getAtoms()) {
        ret += atom_grids[a.getXSType()].getEnergy(a);
      }
      return ret;
    }
  public:
    explicit Optimizer_Grid(const std::vector<fragdock::AtomInterEnergyGrid>& atom_grids) : atom_grids(atom_grids) {}
    fltype optimize(Molecule &mol) const;
  };
}

#endif
