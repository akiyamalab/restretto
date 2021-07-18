#include "common.hpp"
#include "EnergyGrid.hpp"
#include "Molecule.hpp"
#include <vector>
#ifndef ENERGY_CALCULATOR_H_
#define ENERGY_CALCULATOR_H_

namespace fragdock {
  class EnergyCalculator {
  private:
    static const int THRESHOLD = 8;
    static const int PRECI = 1000;
    static const int SZ = THRESHOLD * PRECI;
    static fltype sqr(fltype x) { return x * x; }
    static fltype getEnergy_strict(const Atom &latom, const Atom &ratom);
    std::vector<fltype> pre_calculated_energy;
    bool pre_calculated = false;
    fltype getEnergy(const Atom &latom, const Atom &ratom) const;
  public:
    EnergyCalculator() {}
    EnergyCalculator(fltype rad_scale, fltype threthold);
    fltype getEnergy(const Atom &latom, const Molecule &receptor) const;
    fltype getEnergy(const Molecule &ligand, const Molecule &receptor) const;
    fltype getIntraEnergy(const Molecule &ligand) const;

    static fltype gauss1(const Molecule &ligand, const Molecule &receptor);
    static fltype gauss2(const Molecule &ligand, const Molecule &receptor);
    static fltype repulsion(const Molecule &ligand, const Molecule &receptor);
    static fltype hydrophobic(const Molecule &ligand, const Molecule &receptor);
    static fltype Hydrogen(const Molecule &ligand, const Molecule &receptor);
    static fltype getIntraEnergy_strict(const Molecule &ligand);

  };
}
#endif
