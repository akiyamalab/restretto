#include "common.hpp"
#include "EnergyGrid.hpp"
#include "AtomEnergyGrid.hpp"
#include "Molecule.hpp"
#include <vector>
#ifndef GRADIENT_CALCULATOR_H_
#define GRADIENT_CALCULATOR_H_

namespace fragdock {
  class GradientCalculator {
  private:
    static const int THRESHOLD = 8;
    static const int PRECI = 1000;
    static const int SZ = THRESHOLD * PRECI;
    static fltype sqr(fltype x) { return x * x; }
    static fltype getGradient_strict(const Atom &latom, const Atom &ratom);
    std::vector<fltype> pre_calculated_grad;
    bool pre_calculated = false;
    Vector3d getGradient(const Atom &latom, const Atom &ratom) const;
  public:
    GradientCalculator() {}
    GradientCalculator(fltype rad_scale, fltype threthold);
    Vector3d getGradient(const Atom &latom, const Molecule &receptor) const;
    Vector3d getGradient(const Molecule &ligand, const Molecule &receptor) const;
    static std::vector<std::vector<AtomEnergyGrid> > readGradientGrids(const std::string& grid_folder);
  };
}
#endif
