#include "InterEnergyGrid.hpp"
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <openbabel/griddata.h>
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

  void InterEnergyGrid::parseGrid(const std::string& filename) {
    using namespace std;
    ifstream ifs(filename.c_str(), ios::binary);
    if(!ifs) {
      cerr << "InterEnergyGrid::parseGrid() : file could not open. " << filename << endl;
      abort();
    }
    parseGrid(ifs);
    ifs.close();
  }
  
  void InterEnergyGrid::parseGrid(std::ifstream& ifs) {
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

  void InterEnergyGrid::parseDx(const std::string& filename) {
    OpenBabel::OBConversion conv;
    OpenBabel::OBMol mol;

    if (!conv.SetInFormat("dx")) {
      std::cerr << "InterEnergyGrid::parseDx() : could not set OpenDX format." << std::endl;
      abort();
    }
    if (!conv.ReadFile(&mol, filename)) {
      std::cerr << "InterEnergyGrid::parseDx() : failed to read file " << filename << std::endl;
      abort();
    }

    // Some OpenBabel installations expose only non-template GetData() overloads.
    // Use the typed enum-based lookup and dynamic_cast to obtain OBGridData.
    OpenBabel::OBGenericData* gd = mol.GetData(OpenBabel::OBGenericDataType::GridData);
    OpenBabel::OBGridData* grid = nullptr;
    if (gd) grid = dynamic_cast<OpenBabel::OBGridData*>(gd);
    if (!grid) {
      std::cerr << "InterEnergyGrid::parseDx() : no grid data found in file " << filename << std::endl;
      abort();
    }

    grid->GetNumberOfPoints(num.x, num.y, num.z);

    OpenBabel::vector3 sx, sy, sz;
    grid->GetAxes(sx, sy, sz);
    pitch = Point3d<fltype>(sx.GetX(), sy.GetY(), sz.GetZ());
    
    OpenBabel::vector3 origin = grid->GetOriginVector();
    center = Point3d<fltype>(origin.GetX(), origin.GetY(), origin.GetZ()) + pitch * (num - 1) / 2;

    initInterEnergy();

    for (int x = 0; x < num.x; x++) {
      for (int y = 0; y < num.y; y++) {
        for (int z = 0; z < num.z; z++) {
          fltype val = grid->GetValue(x, y, z);
          setInterEnergy(x, y, z, val);
        }
      }
    }
  }

  void InterEnergyGrid::parse(const std::string& filename) {
    using namespace std;
    auto pos = filename.find_last_of(".");
    if (pos == string::npos) {
      cerr << "InterEnergyGrid::parse() : unknown file extension. " << filename << endl;
      abort();
    }
    string extention = filename.substr(pos);

    if (extention == ".grid") {
      parseGrid(filename);
    } else if (extention == ".dx") {
      parseDx(filename);
    } else {
      cerr << "InterEnergyGrid::parse() : unknown file extension. " << filename << endl;
      abort();
    }
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
      abort();
    }
    writeFile(ofs);
    ofs.close();
  }
}
