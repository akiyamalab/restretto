#include "FragmentEnergyGrid.hpp"

namespace fragdock {
  FragmentEnergyGrid::FragmentEnergyGrid(const Fragment& orig_frag,
                                         const std::vector<Vector3d>& rot_angles,
                                         const std::vector<AtomEnergyGrid>& atom_grids,
                                         const EnergyGrid& distance_grid,
                                         const int top_fgrid_num) : frag_id(orig_frag.getId()),
                                                                            temp_id(orig_frag.gettempind()) {
    using namespace std;
    if(atom_grids.empty()) {
      cerr << "atom_grids is empty" << endl;
      return;
    }
    assert(orig_frag.getCenter().abs() < 1e-4);
    const Point3d<int>& num = atom_grids[0].getNum();
    grid = EnergyGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), num, LIMIT_ENERGY); //atom_grids[0].getPitch()のデフォルト:0.25 (textgrid.3.in)
    
    utils::MinValuesVector<XYZ_param> mvv_good_xyzs = utils::MinValuesVector<XYZ_param>(top_fgrid_num);
    int rotsz = rot_angles.size();
    if (orig_frag.size() <= 1) rotsz = 1;

    fltype radius = orig_frag.getRadius();
    
    // RXYZ -> XYZR の変更に伴い，先に各回転操作を行っておく．
    vector<Fragment> frags(rotsz);
    for(int rot_id = 0; rot_id < rotsz; rot_id++) {
      frags[rot_id] = orig_frag;
      frags[rot_id].rotate(rot_angles[rot_id]);
    }
    
    // RXYZ -> XYZR のループ構造に変更（この方がF2の準備にとって都合が良いため）
    for(int x = 0; x < num.x; x++) {
      for(int y = 0; y < num.y; y++) {
        for(int z = 0; z < num.z; z++) {
          for(int rot_id = 0; rot_id < rotsz; rot_id++) {
            Fragment& frag = frags[rot_id];
            const vector<Atom>& atoms = frag.getAtoms();
            
            if (rot_id & 7) {
              fltype collision = 2;
              fltype far = 6;
              fltype dist = distance_grid.getEnergy(x, y, z); //(x, y, z)と，(x, y, z)に1番近いタンパク質の原子(水素原子以外であれば何でもOK)との距離
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
          
          mvv_good_xyzs.push(XYZ_param(x, y, z, grid.getEnergy(x, y, z)));
          
        }
      }
    }
    
    good_xyzs = mvv_good_xyzs.getValues();
    sort(good_xyzs.begin(), good_xyzs.end());
  }


  FragmentEnergyGrid::FragmentEnergyGrid(const int frag_id,
                                         const std::string& filename,
                                         const int top_fgrid_num) :
                                         frag_id(frag_id) {
    using namespace std;
    parse(filename);

    utils::MinValuesVector<XYZ_param> mvv_good_xyzs = utils::MinValuesVector<XYZ_param>(top_fgrid_num);

    const Point3d<int>& num = grid.getNum();
    for(int x = 0; x < num.x; x++)
    for(int y = 0; y < num.y; y++)
    for(int z = 0; z < num.z; z++)
      mvv_good_xyzs.push(XYZ_param(x, y, z, grid.getEnergy(x, y, z)));


    good_xyzs = mvv_good_xyzs.getValues();
    sort(good_xyzs.begin(), good_xyzs.end());
  }

  void FragmentEnergyGrid::parse(const std::string& filename) {
    using namespace std;
    ifstream ifs(filename.c_str(), std::ios::binary);
    if(!ifs) {
      cerr << "FragmentEnergyGrid::parse() : file could not open. " << filename << endl;
      return;
    }
    grid.parse(ifs);
    // grids.resize(rot_size);
    // for (int i = 0; i < rot_size; i++)
    //  grids[i].parse(ifs);

    ifs.close();
  }

  void FragmentEnergyGrid::writeFile(const std::string& filename) const {
    using namespace std;
    ofstream ofs(filename.c_str(), std::ios::binary);
    if(!ofs){
      cerr << "EnergyGrid::WriteFile() : file could not open. " << filename << endl;
      return ;
    }
    grid.writeFile(ofs);
    // for(int i = 0; i < grids.size(); i++)
    //   grids[i].writeFile(ofs);

    ofs.close();
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
