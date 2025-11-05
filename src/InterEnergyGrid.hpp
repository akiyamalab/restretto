#include "common.hpp"
#include "Point3d.hpp"
#include "Vector3d.hpp"
#include "utils.hpp"
#include <fstream>
#include <iostream>

#ifndef _ENERGY_GRID_H_
#define _ENERGY_GRID_H_

namespace fragdock {

  /**
   * A 3D grid of energy values.
  */
  class InterEnergyGrid {
  private:
    /* initialize the energy grid with a real value ini */
    void initInterEnergy(fltype ini) { this->grid = std::vector<fltype>(num.x*num.y*num.z, ini); }
    void initInterEnergy() { initInterEnergy(INF_ENERGY); }

  protected:
    int max(int x, int y) const { return x>y?x:y; }
    int min(int x, int y) const { return x<y?x:y; }
    std::vector<fltype> grid; /* A raw 3D grid of energy values */
    Point3d<fltype> center; /* The center coordinate of a grid */
    Point3d<fltype> pitch; /* The grid pitch for each dimension */
    Point3d<int> num; /* The numbers of grid points for each dimension */

  public:
    InterEnergyGrid() {}
    InterEnergyGrid(const Point3d<fltype>& center, const Point3d<fltype>& pitch, const Point3d<int>& num)
      : center(center), pitch(pitch), num(num) { initInterEnergy(); }
    InterEnergyGrid(const Point3d<fltype>& center, const Point3d<fltype>& pitch, const Point3d<int>& num, fltype ini)
      : center(center), pitch(pitch), num(num) { initInterEnergy(ini); }
    ~InterEnergyGrid() {}

    /* set an energy for grid indice [x_ind,y_ind,z_ind]. Note that the x,y,z are not coordinates */
    void setInterEnergy(int x_ind, int y_ind, int z_ind, fltype val) { grid[(x_ind*num.y+y_ind)*num.z+z_ind] = val; }

    /* add an energy for grid indice [x_ind,y_ind,z_ind]. Note that the x,y,z are not coordinates */
    void addEnergy(int x_ind, int y_ind, int z_ind, fltype val) { grid[(x_ind*num.y+y_ind)*num.z+z_ind] += val; }


    /* get an energy of grid indice [x_ind,y_ind,z_ind]. Note that the x,y,z are not coordinates*/
    fltype getInterEnergy(int x_ind, int y_ind, int z_ind) const;

    /* get an energy of coordinate pos. */
    fltype getInterEnergy(const Vector3d &pos) const;

    /* get grid pitches for x, y, and z dimensions */
    const Point3d<fltype>& getPitch() const { return pitch; }

    /* get center position of the grid */
    const Point3d<fltype>& getCenter() const { return center; }

    /* get the numbers of grid points */
    const Point3d<int>& getNum() const { return num; }

    /* convert an x coordinate into x_ind of this grid */
    int convertX(const Vector3d &vec) const { return utils::round( (vec.x-center.x)/pitch.x + (num.x-1)/2); }

    /* convert a y coordinate into y_ind of this grid */
    int convertY(const Vector3d &vec) const { return utils::round( (vec.y-center.y)/pitch.y + (num.y-1)/2); }

    /* convert a z coordinate into z_ind of this grid */
    int convertZ(const Vector3d &vec) const { return utils::round( (vec.z-center.z)/pitch.z + (num.z-1)/2); }

    /* convert indices x_ind, y_ind, and z_ind to a 3d coordinate */
    Vector3d convert(int x_ind, int y_ind, int z_ind) const { 
      return Vector3d((x_ind - (num.x-1)/2)*pitch.x + center.x,
								      (y_ind - (num.y-1)/2)*pitch.y + center.y,
								      (z_ind - (num.z-1)/2)*pitch.z + center.z); 
    }
    /* convert grid point indices to a 3d coordinate */
    Vector3d convert(const Point3d<int>& ind) const { return convert(ind.x, ind.y, ind.z); }

    /* Parse an energy grid data */
    void parseGrid(std::ifstream& ifs);

    /* Parse an energy OpenDX data */
    void parseDx(const std::string& filename);

    /* Parse an energy grid file */
    void parse(const std::string& filename);

    /* Write this energy grid data to a stream */
    void writeFile(std::ofstream& ofs) const;

    /* Write this energy grid data to a file */
    void writeFile(const std::string& filename) const;
  };
}

#endif