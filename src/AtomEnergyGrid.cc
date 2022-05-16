#include "AtomEnergyGrid.hpp"

namespace fragdock {
  std::vector<AtomEnergyGrid> AtomEnergyGrid::readAtomGrids(const std::string& grid_folder){
    std::vector<AtomEnergyGrid> energy_grids(XS_TYPE_SIZE);

    for(int i = 0; i < XS_TYPE_SIZE; i++){
      std::string grid_filename = grid_folder + "/" + xs_strings[i] + ".grid";
      energy_grids[i] = AtomEnergyGrid(grid_filename, i);
    }
    return energy_grids;
  }
  std::vector<AtomEnergyGrid> AtomEnergyGrid::makeAtomGrids(const Point3d<fltype>& center,
                                                            const Point3d<fltype>& pitch,
                                                            const Point3d<int>& num,
                                                            const Molecule& receptor_mol,
                                                            const EnergyCalculator& ec) {

    std::vector<AtomEnergyGrid> energy_grids(XS_TYPE_SIZE);

    for(int i = 0; i < XS_TYPE_SIZE; i++) {
      energy_grids[i] = AtomEnergyGrid(center, pitch, num, i);
      const Vector3d temp(0,0,0);
      Atom atom(0, temp, i);

      for(int x=0; x<energy_grids[i].getNum().x; x++) {
        for(int y=0; y<energy_grids[i].getNum().y; y++) {
          for(int z=0; z<energy_grids[i].getNum().z; z++) {
            // fragdock::Atom a = atom;
            // fragdock::Molecule mol(std::vector<fragdock::Atom>(1, atom));
            atom.setPos(energy_grids[i].convert(x, y, z)); //atomを実際の座標へ移動
            energy_grids[i].setEnergy(x, y, z, ec.getEnergy(atom, receptor_mol));
          }
        }
      }

    }
    return energy_grids;
  }
}
