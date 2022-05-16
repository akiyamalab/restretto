#include "EnergyGrid.hpp"
// #define TRILINEAR

namespace fragdock {
  fltype EnergyGrid::getEnergy(int x, int y, int z) const { //(配列上における)インデックスが与えられた場合(必ず格子点)．フラグメントグリッドのアクセス方法．
    if (x < 0 or y < 0 or z < 0 or x >= num.x or y >= num.y or z >= num.z) return LIMIT_ENERGY;
    // x = max(0, min(x, num.x - 1));
    // y = max(0, min(y, num.y - 1));
    // z = max(0, min(z, num.z - 1));
    return grid[(x*num.y+y)*num.z+z];
  }

  fltype EnergyGrid::getEnergy(const Vector3d &pos) const { //(配列上における)インデックスではなく具体的な座標が与えられた場合(格子点上とは限らない)．原子グリッドのアクセス方法(ただデフォルトでは，54行目にあるように近傍の座標の配列上におけるインデックスを実引数として5行目を呼んでいるので，結局はフラグメントグリッドと値の取得の仕方は同じ)．
    fltype energy = 0.0;

#ifdef TRILINEAR
    Point3d<fltype> real((pos.x-center.x)/pitch.x + (num.x-1)/2,
			 (pos.y-center.y)/pitch.y + (num.y-1)/2,
			 (pos.z-center.z)/pitch.z + (num.z-1)/2);

    //trilinear interpolation トリリニア補間．単に近くの格子点に近似したりするといった近似よりもしっかりとした近似(IPSJ 2015)
    int x = floor(real.x);
    if(x < 0) x = 0;
    else if(x >= num.x-1) x = num.x-2;
    fltype frac_x = real.x - x;
    frac_x = (frac_x<0 ? 0 : (frac_x>1?1:frac_x) );
    fltype revf_x = 1 - frac_x;

    int y = floor(real.y);
    if(y < 0) y = 0;
    else if(y >= num.y-1)   y = num.y-2;
    fltype frac_y = real.y - y;
    frac_y = (frac_y<0 ? 0 : (frac_y>1?1:frac_y) );
    fltype revf_y = 1 - frac_y;

    int z = floor(real.z);
    if(z < 0) z = 0;
    else if(z >= num.z-1)   z = num.z-2;
    fltype frac_z = real.z - z;
    frac_z = (frac_z<0 ? 0 : (frac_z>1?1:frac_z) );
    fltype revf_z = 1 - frac_z;


    energy += getEnergy(x + 1, y + 1, z + 1) * frac_x * frac_y * frac_z;
    energy += getEnergy(x + 1, y + 1, z    ) * frac_x * frac_y * revf_z;
    energy += getEnergy(x + 1, y    , z + 1) * frac_x * revf_y * frac_z;
    energy += getEnergy(x + 1, y    , z    ) * frac_x * revf_y * revf_z;
    energy += getEnergy(x    , y + 1, z + 1) * revf_x * frac_y * frac_z;
    energy += getEnergy(x    , y + 1, z    ) * revf_x * frac_y * revf_z;
    energy += getEnergy(x    , y    , z + 1) * revf_x * revf_y * frac_z;
    energy += getEnergy(x    , y    , z    ) * revf_x * revf_y * revf_z;
#endif
#ifndef TRILINEAR //デフォルトはこちら．トリリニア補間でない場合．if"n"def に注意．
    energy = getEnergy(convertX(pos), convertY(pos), convertZ(pos)); //これは単に近くの格子点一つに近似してそこのスコアを返す．
#endif
    return energy;
  }
  void EnergyGrid::parse(std::ifstream& ifs) {
    ifs.read(reinterpret_cast<char*>(&center), sizeof(center));
    ifs.read(reinterpret_cast<char*>(&pitch), sizeof(pitch));
    ifs.read(reinterpret_cast<char*>(&num), sizeof(num));
    initEnergy();
    for(int x = 0; x < num.x; x++) {
      for(int y = 0; y < num.y; y++) {
        for(int z = 0; z < num.z; z++) {
          fltype val;
          ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
          setEnergy(x, y, z, val);
        }
      }
    }
  }

  void EnergyGrid::parse(const std::string& filename) {
    using namespace std;
    ifstream ifs(filename.c_str(), ios::binary);
    if(!ifs) {
      cerr << "EnergyGrid::parse() : file could not open. " << filename << endl;
      return;
    }
    parse(ifs);
    ifs.close();
  }
  void EnergyGrid::writeFile(std::ofstream& ofs) const {
    ofs.write(reinterpret_cast<const char*>(&center), sizeof(center));
    ofs.write(reinterpret_cast<const char*>(&pitch),  sizeof(pitch) );
    ofs.write(reinterpret_cast<const char*>(&num),    sizeof(num)   );

    for(int x = 0; x < num.x; x++) {
      for(int y = 0; y < num.y; y++) {
        for(int z = 0; z < num.z; z++) {
          fltype val = getEnergy(x, y, z);
          ofs.write(reinterpret_cast<const char*>(&val), sizeof(val));
        }
      }
    }
  }
  void EnergyGrid::writeFile(const std::string& filename) const {
    using namespace std;
    ofstream ofs(filename.c_str(), ios::binary);
    if(!ofs){
      cerr << "EnergyGrid::WriteFile() : file could not open. " << filename << endl;
      return ;
    }
    writeFile(ofs);
    ofs.close();
  }
}
