#include "common.hpp"
#include "EnergyGrid.hpp"
#include "Molecule.hpp"
#include <vector>
#ifndef ENERGY_CALCULATOR_H_
#define ENERGY_CALCULATOR_H_

namespace fragdock {
  /**
   * Class for pre-calculating interaction energy between two atoms.
   * It is used for rapid calculation of energy with slight loss of accuracy.
  */
  class EnergyCalculator {
  private:
    /* calculate energy without using pre-calculation results*/
    static fltype getEnergy_strict(const Atom &latom, const Atom &ratom); 

    std::vector<fltype> pre_calculated_energy;
    bool pre_calculated = false;

    /* look-up an interaction energy between two atoms */
    fltype getEnergy(const Atom &latom, const Atom &ratom) const; 

    /**
     * gauss_1 term of AutoDock Vina
     * @param t1 type of atom 1
     * @param t2 type of atom 2
     * @param d surface distance between atom 1 and atom 2 
    */
    static fltype gauss1(int t1, int t2, fltype d);
    static fltype gauss2(int t1, int t2, fltype d);
    static fltype repulsion(int t1, int t2, fltype d);
    static fltype hydrophobic(int t1, int t2, fltype d);
    static fltype hydrogenBond(int t1, int t2, fltype d);
    const std::vector<fltype> term_weights = {
      -0.035579, // gauss1
      -0.005156, // gauss2
      0.840245,  // repulsion
      -0.035069, // hydrophobic
      -0.587439  // hydrogen
    };
  public:
    EnergyCalculator() {}

    /**
     * @param rad_scale scale factor of radius of atoms
    */
    EnergyCalculator(fltype rad_scale);

    /* calculate an interaction energy between an ligand atom and all protein atoms */
    fltype getEnergy(const Atom &latom, const Molecule &receptor) const; 

    /* calculate an interaction energy between all ligand atoms and all protein atoms */
    fltype getEnergy(const Molecule &ligand, const Molecule &receptor) const;

    /* calculate an intra energy of a ligand */
    fltype calcIntraEnergy(const Molecule &ligand) const;

    static fltype gauss1(const Molecule &ligand, const Molecule &receptor);
    static fltype gauss2(const Molecule &ligand, const Molecule &receptor);
    static fltype repulsion(const Molecule &ligand, const Molecule &receptor);
    static fltype hydrophobic(const Molecule &ligand, const Molecule &receptor);
    static fltype hydrogenBond(const Molecule &ligand, const Molecule &receptor);
    static fltype getIntraEnergy_strict(const Molecule &ligand);

  };
}
#endif
