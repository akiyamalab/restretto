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

  void InterEnergyGrid::parseDx(std::ifstream& ifs) {
    Point3d<fltype> origin;
    pitch = Point3d<fltype>(0.0, 0.0, 0.0);

    std::string line;
    while (std::getline(ifs, line)) {
      std::istringstream iss(line);
      std::string token;
      iss >> token;

      if (token == "object") {
        iss >> token;
        if (token == "1") {
          std::string tmp;
          iss >> tmp >> tmp >> tmp;
          iss >> num.x >> num.y >> num.z;
        } else if (token == "3") {
          break;
        }
      } else if (token == "origin") {
        fltype x, y, z;
        iss >> x >> y >> z;
        origin = Point3d<fltype>(x, y, z);
      } else if (token == "delta") {
        fltype x, y, z;
        iss >> x >> y >> z;
        pitch.x = std::max(pitch.x, x);
        pitch.y = std::max(pitch.y, y);
        pitch.z = std::max(pitch.z, z);
      }
    }
    center = origin + pitch * (num-1) / 2;
    initInterEnergy();
    for (int x = 0; x < num.x; x++) {
      for (int y = 0; y < num.y; y++) {
        for (int z = 0; z < num.z; z++) {
          fltype val;
          ifs >> val;
          setInterEnergy(x, y, z, val);
        }
      }
    }
  }

  void InterEnergyGrid::parse(const std::string& filename) {
    using namespace std;
    ifstream ifs;
    auto pos = filename.find_last_of(".");
    if (pos == string::npos) {
      cerr << "InterEnergyGrid::parse() : unknown file extension. " << filename << endl;
      return;
    }
    string extention = filename.substr(pos);

    if (extention == ".grid") {
      ifs.open(filename.c_str(), ios::binary);
      if(!ifs) {
        cerr << "InterEnergyGrid::parse() : file could not open. " << filename << endl;
        return;
      }
      parseGrid(ifs);
      ifs.close();
    } else if (extention == ".dx") {
      ifs.open(filename.c_str());
      if(!ifs) {
        cerr << "InterEnergyGrid::parse() : file could not open. " << filename << endl;
        return;
      }
      parseDx(ifs);
      ifs.close();
    } else {
      cerr << "InterEnergyGrid::parse() : unknown file extension. " << filename << endl;
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
      return ;
    }
    writeFile(ofs);
    ofs.close();
  }
}
