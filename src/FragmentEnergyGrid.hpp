#include "common.hpp"
#include "InterEnergyGrid.hpp"
// #include "EnergyCalculator.hpp"
#include "AtomEnergyGrid.hpp"
#include "Fragment.hpp"
#include "Point3d.hpp"

#ifndef FRAGMENT_ENERGY_GRID_H_
#define FRAGMENT_ENERGY_GRID_H_
namespace fragdock {
  class FragmentEnergyGrid /* : public InterEnergyGrid */ {
    InterEnergyGrid grid;
    // void parse(const std::string& filename, int rot_size);
  public:
    int frag_idx;
    FragmentEnergyGrid() { frag_idx = -1; }
    // FragmentEnergyGrid(int frag_id, const std::string& filename, int rot_size)
    // : frag_id(frag_id) { parse(filename, rot_size); }
    FragmentEnergyGrid(const Fragment& orig_frag,
                       const std::vector<Vector3d>& rot_angles,
                       const std::vector<AtomEnergyGrid>& atom_grids,
                       const InterEnergyGrid& distance_grid);
    ~FragmentEnergyGrid() {}
    const InterEnergyGrid& getGrid() const { return grid; }
    // void writeFile(const std::string& filename) const;
  };
}
#endif

