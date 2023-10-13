#include "InterEnergyGrid.hpp"
// #define TRILINEAR

namespace fragdock {
  fltype InterEnergyGrid::getInterEnergy(int x, int y, int z) const {
    if (x < 0 or y < 0 or z < 0 or x >= num.x or y >= num.y or z >= num.z) return LIMIT_ENERGY;
    // x = max(0, min(x, num.x - 1));
    // y = max(0, min(y, num.y - 1));
    // z = max(0, min(z, num.z - 1));
    return grid[(x*num.y+y)*num.z+z];
  }

  fltype InterEnergyGrid::getInterEnergy(const Vector3d &pos) const {
    return getInterEnergy(convertX(pos), convertY(pos), convertZ(pos)); // nearest grid point
  }
  void InterEnergyGrid::parse(std::ifstream& ifs) {
    ifs.read(reinterpret_cast<char*>(&center), sizeof(center));
    ifs.read(reinterpret_cast<char*>(&pitch), sizeof(pitch));
    ifs.read(reinterpret_cast<char*>(&num), sizeof(num));
    initInterEnergy();
    for(int x = 0; x < num.x; x++) {
      for(int y = 0; y < num.y; y++) {
        for(int z = 0; z < num.z; z++) {
          fltype val;
          ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
          setInterEnergy(x, y, z, val);
        }
      }
    }
  }

  void InterEnergyGrid::parse(const std::string& filename) {
    using namespace std;
    ifstream ifs(filename.c_str(), ios::binary);
    if(!ifs) {
      cerr << "InterEnergyGrid::parse() : file could not open. " << filename << endl;
      return;
    }
    parse(ifs);
    ifs.close();
  }
  void InterEnergyGrid::writeFile(std::ofstream& ofs) const {
    ofs.write(reinterpret_cast<const char*>(&center), sizeof(center));
    ofs.write(reinterpret_cast<const char*>(&pitch),  sizeof(pitch) );
    ofs.write(reinterpret_cast<const char*>(&num),    sizeof(num)   );

    for(int x = 0; x < num.x; x++) {
      for(int y = 0; y < num.y; y++) {
        for(int z = 0; z < num.z; z++) {
          fltype val = getInterEnergy(x, y, z);
          ofs.write(reinterpret_cast<const char*>(&val), sizeof(val));
        }
      }
    }
  }
  void InterEnergyGrid::writeFile(const std::string& filename) const {
    using namespace std;
    ofstream ofs(filename.c_str(), ios::binary);
    if(!ofs){
      cerr << "InterEnergyGrid::WriteFile() : file could not open. " << filename << endl;
      return ;
    }
    writeFile(ofs);
    ofs.close();
  }
}
