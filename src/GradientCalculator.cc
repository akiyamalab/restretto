#include "GradientCalculator.hpp"

namespace fragdock {
  GradientCalculator::GradientCalculator(fltype rad_scale, fltype threthold) {
    pre_calculated_grad = std::vector<fltype>(XS_TYPE_SIZE * XS_TYPE_SIZE * SZ);
    for (int t1 = 0; t1 < XS_TYPE_SIZE; ++t1) {
      for (int t2 = 0; t2 < XS_TYPE_SIZE; ++t2) {
        fltype rad = (xs_radius(t1) + xs_radius(t2)) * rad_scale;
        for (int i = 0; i < SZ; ++i) {
          fltype r = i / (fltype)PRECI;
          fltype d = r - rad;

          pre_calculated_grad[(t1 * XS_TYPE_SIZE + t2) * SZ + i]
              = (-0.035579) * (-8 * d * exp(-sqr((2 * d)))) * (t1 != XS_TYPE_DUMMY && t2 != XS_TYPE_DUMMY ? 1 : 0)
              + (-0.005156) * (-(d - 3.0) * 0.5)*exp(-sqr((d - 3.0) * 0.5)) * (t1 != XS_TYPE_DUMMY && t2 != XS_TYPE_DUMMY ? 1 : 0)
              + ( 0.840245) * (d > threthold ? 0.0 : 2 * d)
              // + ( 0.840245) * (d > threthold ? 0.0 : sqr(d - threthold))
              // + ( 0.840245) * (d > 0.0 ? 0.0 : ((threthold < -1e-3 && d > threthold) ? sqr(d * d / threthold) : sqr(d)))
              + (-0.035069) * ((xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 0.0 : - 1.0)) : 0.0)
              + (-0.587439) * ((xs_hbond(t1, t2)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 0.0 : -1.428571)): 0.0);
        }
      }
    }
    pre_calculated = true;
  }

  Vector3d GradientCalculator::getGradient(const Atom &latom, const Molecule &receptor) const {
    Vector3d sum_grad(0.0, 0.0, 0.0);
    for(int i = 0; i < receptor.size(); i++) {
      const Atom &ratom = receptor.getAtom(i);
      if (ratom.getXSType() != XS_TYPE_H)
        sum_grad += getGradient(latom, ratom);
    }
    return sum_grad;
  }

  Vector3d GradientCalculator::getGradient(const Atom &latom, const Atom &ratom) const {
    assert(pre_calculated);
    int t1 = latom.getXSType();
    int t2 = ratom.getXSType();
    fltype r = (latom - ratom).abs();
    if (r > THRESHOLD) return Vector3d(0.0, 0.0, 0.0);

    int idx = (int)(r * PRECI);
    if (idx >= SZ) return Vector3d(0.0, 0.0, 0.0);

    fltype derE_derr = pre_calculated_grad[(t1 * XS_TYPE_SIZE + t2) * SZ + idx];
    assert(r >= 0);
    Vector3d grad = (r > 0) ? (latom - ratom) / r * derE_derr 
                            : Vector3d(INF_ENERGY, INF_ENERGY, INF_ENERGY);
    return grad;
  }

  std::vector<std::vector<AtomEnergyGrid> > GradientCalculator::readGradientGrids(const std::string& grid_folder){
    std::vector<std::vector<AtomEnergyGrid> > gradient_grids(XS_TYPE_SIZE, std::vector<AtomEnergyGrid>(3));

    for(int i = 0; i < XS_TYPE_SIZE; i++){
      std::string grid_filename_x = grid_folder + "/" + xs_strings[i] + "_gradient_x.grid";
      std::string grid_filename_y = grid_folder + "/" + xs_strings[i] + "_gradient_y.grid";
      std::string grid_filename_z = grid_folder + "/" + xs_strings[i] + "_gradient_z.grid";
      
      gradient_grids[i][0] = AtomEnergyGrid(grid_filename_x, i);
      gradient_grids[i][1] = AtomEnergyGrid(grid_filename_y, i);
      gradient_grids[i][2] = AtomEnergyGrid(grid_filename_z, i);
    }
    return gradient_grids;
  }
}
