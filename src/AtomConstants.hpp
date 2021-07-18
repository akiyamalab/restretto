#include "common.hpp"

#include <string>

#ifndef ATOM_CONSTANTS_H_
#define ATOM_CONSTANTS_H_


// X-Score
const int XS_TYPE_C_H   =  0;
const int XS_TYPE_C_P   =  1;
const int XS_TYPE_N_P   =  2;
const int XS_TYPE_N_D   =  3;
const int XS_TYPE_N_A   =  4;
const int XS_TYPE_N_DA  =  5;
const int XS_TYPE_O_P   =  6;
const int XS_TYPE_O_D   =  7;
const int XS_TYPE_O_A   =  8;
const int XS_TYPE_O_DA  =  9;
const int XS_TYPE_S_P   = 10;
const int XS_TYPE_P_P   = 11;
const int XS_TYPE_F_H   = 12;
const int XS_TYPE_Cl_H  = 13;
const int XS_TYPE_Br_H  = 14;
const int XS_TYPE_I_H   = 15;
const int XS_TYPE_Met_D = 16;
const int XS_TYPE_OTHER = 17;
const int XS_TYPE_DUMMY = 18;

const int XS_TYPE_SIZE  = 19;

const int XS_TYPE_H     = 20;

const std::string xs_strings[] = {
  "C_H",
  "C_P",
  "N_P",
  "N_D",
  "N_A",
  "N_DA",
  "O_P",
  "O_D",
  "O_A",
  "O_DA",
  "S_P",
  "P_P",
  "F_H",
  "Cl_H",
  "Br_H",
  "I_H",
  "Met_D",
  "Other",
  "Dummy"
};


const fltype xs_vdw_radii[] = {
  1.9, // C_H
  1.9, // C_P
  1.8, // N_P
  1.8, // N_D
  1.8, // N_A
  1.8, // N_DA
  1.7, // O_P
  1.7, // O_D
  1.7, // O_A
  1.7, // O_DA
  2.0, // S_P
  2.1, // P_P
  1.5, // F_H
  1.8, // Cl_H
  2.0, // Br_H
  2.2, // I_H
  1.2, // Met_D
  2.0, // Other
  1.5  // Dummy
};

// inline fltype xs_radius(int xs, fltype scale) {
//   assert(xs < XS_TYPE_SIZE);
//   return xs_vdw_radii[xs] * scale;
// }
inline fltype xs_radius(int xs) {
  assert(xs < XS_TYPE_SIZE);
  return xs_vdw_radii[xs];
}

inline std::string xs_name(int xs) {
  if (xs == XS_TYPE_H) return "H";
  // if (xs == XS_TYPE_Dummy) return "Dummy";
  assert(xs < XS_TYPE_SIZE);
  return xs_strings[xs];
}

inline bool xs_is_hydrophobic(int xs) {
  return xs == XS_TYPE_C_H ||
         xs == XS_TYPE_F_H ||
         xs == XS_TYPE_Cl_H ||
         xs == XS_TYPE_Br_H ||
         xs == XS_TYPE_I_H;
}

inline bool xs_is_acceptor(int xs) {
  return xs == XS_TYPE_N_A ||
         xs == XS_TYPE_N_DA ||
         xs == XS_TYPE_O_A ||
         xs == XS_TYPE_O_DA;
}

inline bool xs_is_donor(int xs) {
  return xs == XS_TYPE_N_D ||
         xs == XS_TYPE_N_DA ||
         xs == XS_TYPE_O_D ||
         xs == XS_TYPE_O_DA ||
         xs == XS_TYPE_Met_D;
}

inline bool xs_donor_acceptor(int t1, int t2) {
  return xs_is_donor(t1) && xs_is_acceptor(t2);
}

inline bool xs_hbond(int t1, int t2) {
  return xs_donor_acceptor(t1, t2) || xs_donor_acceptor(t2, t1);
}
inline bool xs_is_heavy(int xs) {
  return xs != XS_TYPE_H && xs != XS_TYPE_DUMMY;
}


#endif
