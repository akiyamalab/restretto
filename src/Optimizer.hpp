#include "common.hpp"
#include "AtomInterEnergyGrid.hpp"
#include "Molecule.hpp"
#include "utils.hpp"
#ifndef OPTIMIZER_H_
#define OPTIMIZER_H_
namespace fragdock {
  class Optimizer {
    const Molecule& receptor;
    const fltype max_rmsd;
  public:
    explicit Optimizer(const Molecule& receptor, fltype max_rmsd = 1e10) : receptor(receptor), max_rmsd(max_rmsd) {}
    fltype optimize(Molecule &mol, const EnergyCalculator& ec) const;
  };

  class Optimizer_Grid {
    const std::vector<fragdock::AtomInterEnergyGrid>& atom_grids;
    // const Molecule& receptor;
    const fltype max_rmsd;
    fltype calcInterEnergy(const Molecule &mol) const {
      fltype ret = 0.0;
      for (auto& a : mol.getAtoms()) {
        ret += atom_grids[a.getXSType()].getInterEnergy(a);
      }
      return ret;
    }
  public:
    explicit Optimizer_Grid(const std::vector<fragdock::AtomInterEnergyGrid>& atom_grids, fltype max_rmsd = 1e10) : atom_grids(atom_grids), max_rmsd(max_rmsd) {}
    fltype calcTotalEnergy(const Molecule &mol) const {return calcInterEnergy(mol) + mol.getIntraEnergy();}
    fltype optimize(Molecule &mol) const;
  };
}

#endif
