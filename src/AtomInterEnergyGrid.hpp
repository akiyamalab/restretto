#include "common.hpp"
#include "Atom.hpp"
#include "InterEnergyGrid.hpp"
#include "EnergyCalculator.hpp"

#ifndef ATOM_ENERGY_GRID_H_
#define ATOM_ENERGY_GRID_H_
namespace fragdock {
  /**
   * An energy grid for a single atom type.
  */
  class AtomInterEnergyGrid : public InterEnergyGrid {
    int xs_type;
  public:
    AtomInterEnergyGrid() {}
    AtomInterEnergyGrid(const std::string &filename, int xs_type)
      : xs_type(xs_type) { parse(filename); }
    AtomInterEnergyGrid(const Point3d<fltype>& center, const Point3d<fltype>& pitch, const Point3d<int>& num, int xs_type)
      : InterEnergyGrid(center, pitch, num), xs_type(xs_type) {}
    ~AtomInterEnergyGrid() {}
    int getXSType() const { return xs_type; }
    static std::vector<AtomInterEnergyGrid> readAtomGrids(const std::string& grid_folder);
    static std::vector<std::pair<bool, AtomInterEnergyGrid> > readDxAtomGrids(const std::string& dx_folder);
    static std::vector<AtomInterEnergyGrid> makeAtomGrids(const Point3d<fltype>& center,
                                                     const Point3d<fltype>& pitch,
                                                     const Point3d<int>& num,
                                                     const Molecule& receptor_mol,
                                                     const EnergyCalculator& ec);
  };
}
#endif

