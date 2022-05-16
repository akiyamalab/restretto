#include "common.hpp"
#include "EnergyGrid.hpp"
// #include "EnergyCalculator.hpp"
#include "AtomEnergyGrid.hpp"
#include "Fragment.hpp"

#include "MinValuesVector.hpp"

#ifndef FRAGMENT_ENERGY_GRID_H_
#define FRAGMENT_ENERGY_GRID_H_
namespace fragdock {
  struct XYZ_param {
    Vector3d pos;
    fltype score;
    XYZ_param() : score(INF_ENERGY) {
      pos.x = -1; pos.y = -1; pos.z = -1;
    }
    XYZ_param(int x, int y, int z, fltype score) : score(score) {
      pos.x = x; pos.y = y; pos.z = z;
    }
    bool operator<(const XYZ_param& o) const { return score < o.score; }
    bool operator>(const XYZ_param& o) const { return score > o.score; }
  };
  
  class FragmentEnergyGrid /* : public EnergyGrid */ {
    EnergyGrid grid;
    std::vector<XYZ_param> good_xyzs;
    // void parse(const std::string& filename, int rot_size);
    void parse(const std::string& filename);
  public:
    int frag_id;
    int temp_id;
    FragmentEnergyGrid() { frag_id = temp_id = -1; }
    // FragmentEnergyGrid(int frag_id, const std::string& filename, int rot_size)
    // : frag_id(frag_id) { parse(filename, rot_size); }
    FragmentEnergyGrid(const int frag_id,
                       const std::string& filename,
                       const int top_fgrid_num);
    FragmentEnergyGrid(const Fragment& orig_frag,
                       const std::vector<Vector3d>& rot_angles,
                       const std::vector<AtomEnergyGrid>& atom_grids,
                       const EnergyGrid& distance_grid,
                       const int top_fgrid_num);
    ~FragmentEnergyGrid() {}
    const EnergyGrid& getGrid() const { return grid; }
    void writeFile(const std::string& filename) const;
    const fltype& getBestScore() const { return good_xyzs[0].score; }
    const std::vector<XYZ_param>& getGoodXYZs() const { return good_xyzs; }
  };
}
#endif

