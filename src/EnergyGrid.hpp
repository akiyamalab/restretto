#include "common.hpp"
#include "utils.hpp"
#include "Vector3d.hpp"
#include <fstream>
#include <iostream>
#ifndef _ENERGY_GRID_H_
#define _ENERGY_GRID_H_
namespace fragdock {
  template<class T>
  struct Point3d {
    T x, y, z;
    Point3d() : x(T()), y(T()), z(T()) {}
    Point3d(const T &x, const T &y, const T &z) : x(x), y(y), z(z) {}
  };

  class EnergyGrid {
  protected:
    int max(int x, int y) const { return x>y?x:y; }
    int min(int x, int y) const { return x<y?x:y; }
    std::vector<fltype> grid;
    Point3d<fltype> center;
    Point3d<fltype> pitch;
    Point3d<int> num;

  public:
    EnergyGrid() {}
    EnergyGrid(const Point3d<fltype>& center, const Point3d<fltype>& pitch, const Point3d<int>& num)
      : center(center), pitch(pitch), num(num) { initEnergy(); }
    EnergyGrid(const Point3d<fltype>& center, const Point3d<fltype>& pitch, const Point3d<int>& num, fltype ini)
      : center(center), pitch(pitch), num(num) { initEnergy(ini); }
    ~EnergyGrid() {}
    void setEnergy(int x, int y, int z, fltype val) { grid[(x*num.y+y)*num.z+z] = val; }
    // void setEnergy(const Vector3d &pos, fltype val) { setEnergy(convertX(pos), convertY(pos), convertZ(pos), val); }
    void setEnergy(const std::vector<fltype>& grid) { this->grid = grid; }
    void addEnergy(int x, int y, int z, fltype val) { grid[(x*num.y+y)*num.z+z] += val; }
    void initEnergy(fltype ini) { this->grid = std::vector<fltype>(num.x*num.y*num.z, ini); }
    void initEnergy() { initEnergy(INF_ENERGY); }
    fltype getEnergy(int x, int y, int z) const;
    fltype getEnergy(const Vector3d &pos) const;
    const Point3d<fltype>& getPitch() const { return pitch; }
    const Point3d<fltype>& getCenter() const { return center; }
    const Point3d<int>& getNum() const { return num; }
    // int convertX(const Vector3d &vec) const { return max(0, min(num.x-1, utils::round( (vec.x-center.x)/pitch.x + (num.x-1)/2))); }
    // int convertY(const Vector3d &vec) const { return max(0, min(num.y-1, utils::round( (vec.y-center.y)/pitch.y + (num.y-1)/2))); }
    // int convertZ(const Vector3d &vec) const { return max(0, min(num.z-1, utils::round( (vec.z-center.z)/pitch.z + (num.z-1)/2))); }
    int convertX(const Vector3d &vec) const { return utils::round( (vec.x-center.x)/pitch.x + (num.x-1)/2); } //実際の座標vecから(配列上における)インデックスへの変換
    int convertY(const Vector3d &vec) const { return utils::round( (vec.y-center.y)/pitch.y + (num.y-1)/2); }
    int convertZ(const Vector3d &vec) const { return utils::round( (vec.z-center.z)/pitch.z + (num.z-1)/2); }
    Vector3d convert(int x, int y, int z) const { return Vector3d((x - (num.x-1)/2)*pitch.x + center.x,
								  (y - (num.y-1)/2)*pitch.y + center.y,
								  (z - (num.z-1)/2)*pitch.z + center.z); } //(配列上における)インデックス(x, y, z)から実際の座標への変換
    //convertX,Y,Zとconvertは名前がほとんど同じにもかかわらずやっていることが真逆なので注意．convertと同じ処理をx方向のみに関して行うのがconvertXということ"ではない"ので注意．
    void parse(std::ifstream& ifs);
    void parse(const std::string& filename);
    void writeFile(std::ofstream& ofs) const;
    void writeFile(const std::string& filename) const;
  };
}
#endif

