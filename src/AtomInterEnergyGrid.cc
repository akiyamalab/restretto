#include "AtomInterEnergyGrid.hpp"

namespace fragdock {
  std::vector<AtomInterEnergyGrid> AtomInterEnergyGrid::readAtomGrids(const std::string& grid_folder){
    std::vector<AtomInterEnergyGrid> inter_energy_grids(XS_TYPE_SIZE);

    for(int i = 0; i < XS_TYPE_SIZE; i++){
      std::string grid_filename = grid_folder + "/" + xs_strings[i] + ".grid";
      inter_energy_grids[i] = AtomInterEnergyGrid(grid_filename, i);
    }
    return inter_energy_grids;
  }

  std::vector<AtomInterEnergyGrid> AtomInterEnergyGrid::readDxAtomGrids(const std::string& dx_folder){
    std::vector<AtomInterEnergyGrid> inter_energy_grids = {};

    for (int i = 0; i < XS_TYPE_SIZE; i++) {
      std::string dx_filename = dx_folder + "/" + xs_strings[i] + ".dx";
      if (std::ifstream(dx_filename).good()) {
        inter_energy_grids.push_back(AtomInterEnergyGrid(dx_filename, i));
      }
    }
    return inter_energy_grids;
  }

  std::vector<AtomInterEnergyGrid> AtomInterEnergyGrid::makeAtomGrids(const Point3d<fltype>& center,
                                                            const Point3d<fltype>& pitch,
                                                            const Point3d<int>& num,
                                                            const Molecule& receptor_mol,
                                                            const EnergyCalculator& ec) {

    std::vector<AtomInterEnergyGrid> inter_energy_grids(XS_TYPE_SIZE);

    for(int i = 0; i < XS_TYPE_SIZE; i++) {
      inter_energy_grids[i] = AtomInterEnergyGrid(center, pitch, num, i);
      const Vector3d temp(0,0,0);
      Atom atom(0, temp, i);

      for(int x=0; x<inter_energy_grids[i].getNum().x; x++) {
        for(int y=0; y<inter_energy_grids[i].getNum().y; y++) {
          for(int z=0; z<inter_energy_grids[i].getNum().z; z++) {
            // fragdock::Atom a = atom;
            // fragdock::Molecule mol(std::vector<fragdock::Atom>(1, atom));
            atom.setPos(inter_energy_grids[i].convert(x, y, z));
            inter_energy_grids[i].setInterEnergy(x, y, z, ec.getEnergy(atom, receptor_mol));
          }
        }
      }

    }
    return inter_energy_grids;
  }
}
