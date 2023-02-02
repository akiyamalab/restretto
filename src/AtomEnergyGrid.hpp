#include "common.hpp"
#include "Atom.hpp"
#include "EnergyGrid.hpp"
#include "EnergyCalculator.hpp"

#ifndef ATOM_ENERGY_GRID_H_
#define ATOM_ENERGY_GRID_H_
namespace fragdock {
  /**
   * An energy grid for a single atom type.
  */
  class AtomEnergyGrid : public EnergyGrid {
    int xs_type;
  public:
    AtomEnergyGrid() {}
    AtomEnergyGrid(const std::string &filename, int xs_type)
      : xs_type(xs_type) { parse(filename); }
    AtomEnergyGrid(const Point3d<fltype>& center, const Point3d<fltype>& pitch, const Point3d<int>& num, int xs_type)
      : EnergyGrid(center, pitch, num), xs_type(xs_type) {}
    ~AtomEnergyGrid() {}
    int getXSType() const { return xs_type; }
    static std::vector<AtomEnergyGrid> readAtomGrids(const std::string& grid_folder);
    static std::vector<AtomEnergyGrid> makeAtomGrids(const Point3d<fltype>& center,
                                                     const Point3d<fltype>& pitch,
                                                     const Point3d<int>& num,
                                                     const Molecule& receptor_mol,
                                                     const EnergyCalculator& ec);
  };
}
#endif

