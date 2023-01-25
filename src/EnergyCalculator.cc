#include "EnergyCalculator.hpp"

namespace {
  /* distance threshold */
  static const int THRESHOLD = 8; 

  /* precision of pre-calculation*/
  static const int PRECI = 1000;

  /* list size per atom pair*/
  static const int SZ = THRESHOLD * PRECI;
}

namespace fragdock {
  fltype sqr(fltype x) { return x * x; }

  EnergyCalculator::EnergyCalculator(fltype rad_scale) {
    pre_calculated_energy = std::vector<fltype>(XS_TYPE_SIZE * XS_TYPE_SIZE * SZ);
    for (int t1 = 0; t1 < XS_TYPE_SIZE; ++t1) {
      for (int t2 = 0; t2 < XS_TYPE_SIZE; ++t2) {
        fltype rad = (xs_radius(t1) + xs_radius(t2)) * rad_scale;
        for (int i = 0; i < SZ; ++i) {
          fltype r = i / (fltype)PRECI;
          fltype d = r - rad;

          pre_calculated_energy[(t1 * XS_TYPE_SIZE + t2) * SZ + i]
              = (-0.035579) * gauss1(t1, t2, d)
              + (-0.005156) * gauss2(t1, t2, d)
              + ( 0.840245) * repulsion(t1, t2, d)
              + (-0.035069) * hydrophobic(t1, t2, d)
              + (-0.587439) * Hydrogen(t1, t2, d);
        }
      }
    }
    pre_calculated = true;
  }

  fltype EnergyCalculator::getEnergy(const Atom &latom, const Molecule &receptor) const {
    fltype sum_energy = 0.0;
    for(int i = 0; i < receptor.size(); i++) {
      const Atom &ratom = receptor.getAtom(i);
      if (ratom.getXSType() != XS_TYPE_H)
        sum_energy += getEnergy(latom, ratom);
      if (sum_energy >= LIMIT_ENERGY) return LIMIT_ENERGY;
    }
    return sum_energy;
  }
  fltype EnergyCalculator::getEnergy(const Molecule &ligand, const Molecule &receptor) const {
    fltype sum_energy = 0.0;
    for(int i = 0; i < ligand.size(); i++) {
      const Atom &latom = ligand.getAtom(i);
      if (latom.getXSType() != XS_TYPE_H)
        sum_energy += getEnergy(latom, receptor);
      if (sum_energy >= LIMIT_ENERGY) return LIMIT_ENERGY;
    }
    return sum_energy;
  }

  fltype EnergyCalculator::getEnergy(const Atom &latom, const Atom &ratom) const {
    assert(pre_calculated);
    int t1 = latom.getXSType();
    int t2 = ratom.getXSType();
    fltype r = (latom - ratom).abs();
    if (r > THRESHOLD) return 0;

    int idx = (int)(r * PRECI);
    if (idx >= SZ) return 0;

    return pre_calculated_energy[(t1 * XS_TYPE_SIZE + t2) * SZ + idx];
  }

  fltype EnergyCalculator::getIntraEnergy(const Molecule &ligand) const {
    fltype sum_energy = 0.0;
    std::vector<std::vector<int>> dist = ligand.getGraphDistances();
    for (int i = 0; i < ligand.size(); i++) {
      if (ligand.getAtom(i).getXSType() == XS_TYPE_H) continue;
      for (int j = i + 1; j < ligand.size(); j++) {
        if (ligand.getAtom(j).getXSType() == XS_TYPE_H || dist[ligand.getAtom(i).getId()][ligand.getAtom(j).getId()] < 4) continue;

        sum_energy += getEnergy(ligand.getAtom(i), ligand.getAtom(j));
        if (sum_energy >= LIMIT_ENERGY) return LIMIT_ENERGY;
      }
    }
    return sum_energy;
  }



  fltype EnergyCalculator::gauss1(int t1, int t2, fltype d) {
    if (t1 == XS_TYPE_H || t1 == XS_TYPE_DUMMY) return 0;
    if (t2 == XS_TYPE_H || t2 == XS_TYPE_DUMMY) return 0;
    return exp(-sqr(d * 2));
  }

  fltype EnergyCalculator::gauss2(int t1, int t2, fltype d) {
    if (t1 == XS_TYPE_H || t1 == XS_TYPE_DUMMY) return 0;
    if (t2 == XS_TYPE_H || t2 == XS_TYPE_DUMMY) return 0;
    return exp(-sqr((d - 3.0) * 0.5));
  }

  fltype EnergyCalculator::repulsion(int t1, int t2, fltype d) {
    if (t1 == XS_TYPE_H || t1 == XS_TYPE_DUMMY) return 0;
    if (t2 == XS_TYPE_H || t2 == XS_TYPE_DUMMY) return 0;
    return (d > 0 ? 0.0 : sqr(d));
  }

  fltype EnergyCalculator::hydrophobic(int t1, int t2, fltype d) {
    if (t1 == XS_TYPE_H || t1 == XS_TYPE_DUMMY) return 0;
    if (t2 == XS_TYPE_H || t2 == XS_TYPE_DUMMY) return 0;
    if (!xs_is_hydrophobic(t1) || !xs_is_hydrophobic(t2)) return 0;
    return ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d));
  }

  fltype EnergyCalculator::Hydrogen(int t1, int t2, fltype d) {
    if (t1 == XS_TYPE_H || t1 == XS_TYPE_DUMMY) return 0;
    if (t2 == XS_TYPE_H || t2 == XS_TYPE_DUMMY) return 0;
    if (!xs_hbond(t1, t2)) return 0;
    return ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.428571)));
 }

  fltype EnergyCalculator::gauss1(const Molecule &ligand, const Molecule &receptor) {
    fltype sum_energy = 0.0;
    for (const Atom& la : ligand.getAtoms()) {
      if (la.getXSType() == XS_TYPE_H || la.getXSType() == XS_TYPE_DUMMY) continue;
      for (const Atom& ra : receptor.getAtoms()) {
        if (ra.getXSType() == XS_TYPE_H || ra.getXSType() == XS_TYPE_DUMMY) continue;

        int t1 = la.getXSType();
        int t2 = ra.getXSType();
        fltype r = (la - ra).abs();
        if (r > THRESHOLD) continue;
        fltype d = r - xs_radius(t1) - xs_radius(t2);
        sum_energy += gauss1(t1, t2, d);
      }
    }
    return sum_energy;
  }
  fltype EnergyCalculator::gauss2(const Molecule &ligand, const Molecule &receptor) {
    fltype sum_energy = 0.0;
    for (const Atom& la : ligand.getAtoms()) {
      if (la.getXSType() == XS_TYPE_H || la.getXSType() == XS_TYPE_DUMMY) continue;
      for (const Atom& ra : receptor.getAtoms()) {
        if (ra.getXSType() == XS_TYPE_H || ra.getXSType() == XS_TYPE_DUMMY) continue;

        int t1 = la.getXSType();
        int t2 = ra.getXSType();
        fltype r = (la - ra).abs();
        if (r > THRESHOLD) continue;
        fltype d = r - xs_radius(t1) - xs_radius(t2);
        sum_energy += gauss2(t1, t2, d);
      }
    }
    return sum_energy;
  }
  fltype EnergyCalculator::repulsion(const Molecule &ligand, const Molecule &receptor) {
    fltype sum_energy = 0.0;
    for (const Atom& la : ligand.getAtoms()) {
      if (la.getXSType() == XS_TYPE_H || la.getXSType() == XS_TYPE_DUMMY) continue;
      for (const Atom& ra : receptor.getAtoms()) {
        if (ra.getXSType() == XS_TYPE_H || ra.getXSType() == XS_TYPE_DUMMY) continue;

        int t1 = la.getXSType();
        int t2 = ra.getXSType();
        fltype r = (la - ra).abs();
        if (r > THRESHOLD) continue;
        fltype d = r - xs_radius(t1) - xs_radius(t2);
        sum_energy += repulsion(t1, t2, d);
      }
    }
    return sum_energy;
  }
  fltype EnergyCalculator::hydrophobic(const Molecule &ligand, const Molecule &receptor) {
    fltype sum_energy = 0.0;
    for (const Atom& la : ligand.getAtoms()) {
      if (la.getXSType() == XS_TYPE_H || la.getXSType() == XS_TYPE_DUMMY) continue;
      for (const Atom& ra : receptor.getAtoms()) {
        if (ra.getXSType() == XS_TYPE_H || ra.getXSType() == XS_TYPE_DUMMY) continue;

        int t1 = la.getXSType();
        int t2 = ra.getXSType();
        fltype r = (la - ra).abs();
        if (r > THRESHOLD) continue;
        fltype d = r - xs_radius(t1) - xs_radius(t2);
        sum_energy += hydrophobic(t1, t2, d);
      }
    }
    return sum_energy;
  }
  fltype EnergyCalculator::Hydrogen(const Molecule &ligand, const Molecule &receptor) {
    fltype sum_energy = 0.0;
    for (const Atom& la : ligand.getAtoms()) {
      if (la.getXSType() == XS_TYPE_H || la.getXSType() == XS_TYPE_DUMMY) continue;
      for (const Atom& ra : receptor.getAtoms()) {
        if (ra.getXSType() == XS_TYPE_H || ra.getXSType() == XS_TYPE_DUMMY) continue;

        int t1 = la.getXSType();
        int t2 = ra.getXSType();
        fltype r = (la - ra).abs();
        if (r > THRESHOLD) continue;
        fltype d = r - xs_radius(t1) - xs_radius(t2);
        sum_energy += Hydrogen(t1, t2, d);
      }
    }
    return sum_energy;
  }

  fltype EnergyCalculator::getEnergy_strict(const Atom &latom, const Atom &ratom) {
    int t1 = latom.getXSType();
    int t2 = ratom.getXSType();
    fltype r = (latom - ratom).abs();
    if (r > THRESHOLD) return 0.0;
    fltype d = r - xs_radius(t1) - xs_radius(t2);

    return  (-0.035579) * exp(-sqr(d * 2))
          + (-0.005156) * exp(-sqr((d - 3.0) * 0.5))
          + ( 0.840245) * (d > 0.0 ? 0.0 : sqr(d))
          + (-0.035069) * ((xs_is_hydrophobic(t1) && xs_is_hydrophobic(t2)) ? ((d >= 1.5) ? 0.0 : ((d <= 0.5) ? 1.0 : 1.5 - d)) : 0.0)
          + (-0.587439) * ((xs_hbond(t1, t2)) ? ((d >= 0) ? 0.0 : ((d <= -0.7) ? 1 : d * (-1.428571))): 0.0);
  }


  fltype EnergyCalculator::getIntraEnergy_strict(const Molecule &ligand) {
    fltype sum_energy = 0.0;
    std::vector<std::vector<int>> dist = ligand.getGraphDistances();
    for (int i = 0; i < ligand.size(); i++) {
      if (ligand.getAtom(i).getXSType() == XS_TYPE_H || ligand.getAtom(i).getXSType() == XS_TYPE_DUMMY) continue;
      for (int j = i + 1; j < ligand.size(); j++) {
        if (ligand.getAtom(j).getXSType() == XS_TYPE_H || ligand.getAtom(j).getXSType() == XS_TYPE_DUMMY || dist[ligand.getAtom(i).getId()][ligand.getAtom(j).getId()] < 4) continue;

        sum_energy += getEnergy_strict(ligand.getAtom(i), ligand.getAtom(j));
      }
    }
    return sum_energy;
  }



}
