#include "common.hpp"
#include "FragmentEnergyGrid.hpp"

namespace fragdock {
  FragmentEnergyGrid::FragmentEnergyGrid(const Fragment& orig_frag,
                                         const std::vector<Vector3d>& rot_angles,
                                         const std::vector<AtomEnergyGrid>& atom_grids,
                                         const EnergyGrid& distance_grid) : frag_id(orig_frag.getId()),
                                                                            temp_id(orig_frag.gettempind()) {
    using namespace std;
    if(atom_grids.empty()) {
      cerr << "atom_grids is empty" << endl;
      return;
    }
    assert(orig_frag.getCenter().abs() < EPS);
    const Point3d<int>& num = atom_grids[0].getNum();
    grid = EnergyGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), num, LIMIT_ENERGY);
    int rotsz = rot_angles.size();
    if (orig_frag.size() <= 1) rotsz = 1;

    fltype radius = orig_frag.getRadius();

    for(int rot_id = 0; rot_id < rotsz; rot_id++) {
      Fragment frag = orig_frag;
      frag.rotate(rot_angles[rot_id]);
      const vector<Atom>& atoms = frag.getAtoms();
      for(int x = 0; x < num.x; x++) {
        for(int y = 0; y < num.y; y++) {
          for(int z = 0; z < num.z; z++) {
            if (rot_id & 7) {
              fltype collision = 2;
              fltype far = 6;
              fltype dist = distance_grid.getEnergy(x, y, z);
              if (dist < collision) {
                // collision
                continue;
              }
              if (dist > radius + far) {
                // too far
                continue;
              }
            }

            fltype sum = 0.0;
            for (const Atom& atom : atoms) {
              if (atom.getXSType() == XS_TYPE_H) continue;

              const AtomEnergyGrid& agrid = atom_grids[atom.getXSType()];
              // atom += agrid.convert(x, y, z);
              fltype diff = agrid.getEnergy(atom + agrid.convert(x, y, z));
              // atom -= agrid.convert(x, y, z);

              sum += diff;
              if(sum >= LIMIT_ENERGY) {
                break;
              }
            }
            if (sum < grid.getEnergy(x, y, z)) {
              grid.setEnergy(x, y, z, sum);
            }
          }
        }
      }
    }
  }
  // void FragmentEnergyGrid::parse(const std::string& filename, int rot_size) {
  //   using namespace std;
  //   ifstream ifs(filename.c_str(), std::ios::binary);
  //   if(!ifs) {
  //     cerr << "FragmentEnergyGrid::parse() : file could not open. " << filename << endl;
  //     return;
  //   }
  //   grids.resize(rot_size);
  //   for (int i = 0; i < rot_size; i++)
  //     grids[i].parse(ifs);

  //   ifs.close();
  // }
  // void FragmentEnergyGrid::writeFile(const std::string& filename) const {
  //   using namespace std;
  //   ofstream ofs(filename.c_str(), std::ios::binary);
  //   if(!ofs){
  //     cerr << "EnergyGrid::WriteFile() : file could not open. " << filename << endl;
  //     return ;
  //   }
  //   for(int i = 0; i < grids.size(); i++)
  //     grids[i].writeFile(ofs);

  //   ofs.close();
  // }
}
