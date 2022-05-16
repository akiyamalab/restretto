#include "common.hpp"
#include "AtomEnergyGrid.hpp"
#include "Molecule.hpp"
#include "utils.hpp"
#include <Eigen/Dense>
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
    const std::vector<fragdock::AtomEnergyGrid>& atom_grids;
    const std::vector<std::vector<fragdock::AtomEnergyGrid> >& gradient_grids;
    // const Molecule& receptor;
    fltype calcscore(const Molecule &mol) const {
      fltype ret = 0.0;
      for (auto& a : mol.getAtoms()) {
        ret += atom_grids[a.getXSType()].getEnergy(a);
      }
      return ret;
    }
    
  public:
    //explicit Optimizer_Grid(const std::vector<fragdock::AtomEnergyGrid>& atom_grids) : atom_grids(atom_grids) {}
    explicit Optimizer_Grid(const std::vector<fragdock::AtomEnergyGrid>& atom_grids, const std::vector<std::vector<fragdock::AtomEnergyGrid> >& gradient_grids) : atom_grids(atom_grids), gradient_grids(gradient_grids) {}
    fltype optimize(Molecule &mol) const;
    Eigen::VectorXf calcGradient(Molecule &mol, const Eigen::VectorXf& X);
    fltype optimizeBySteepest(Molecule &mol, const int max_iterations, const fltype step_size);
  };
}

#endif
