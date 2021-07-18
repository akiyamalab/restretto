#include "common.hpp"
#include "EnergyGrid.hpp"
// #include "EnergyCalculator.hpp"
#include "AtomEnergyGrid.hpp"
#include "Fragment.hpp"

#ifndef FRAGMENT_ENERGY_GRID_H_
#define FRAGMENT_ENERGY_GRID_H_
namespace fragdock {
  class FragmentEnergyGrid /* : public EnergyGrid */ {
    EnergyGrid grid;
    // void parse(const std::string& filename, int rot_size);
  public:
    int frag_id;
    int temp_id;
    FragmentEnergyGrid() { frag_id = temp_id = -1; }
    // FragmentEnergyGrid(int frag_id, const std::string& filename, int rot_size)
    // : frag_id(frag_id) { parse(filename, rot_size); }
    FragmentEnergyGrid(const Fragment& orig_frag,
                       const std::vector<Vector3d>& rot_angles,
                       const std::vector<AtomEnergyGrid>& atom_grids,
                       const EnergyGrid& distance_grid);
    ~FragmentEnergyGrid() {}
    const EnergyGrid& getGrid() const { return grid; }
    // void writeFile(const std::string& filename) const;
  };
}
#endif

