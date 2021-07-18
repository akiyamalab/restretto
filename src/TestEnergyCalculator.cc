#include <boost/test/unit_test.hpp>
#include "EnergyCalculator.hpp"

#include <cmath>
#include <vector>

#include "Atom.hpp"
#include "common.hpp"
#include "Molecule.hpp"

BOOST_AUTO_TEST_CASE(CH_CH_interaction){
  using namespace fragdock;
  EnergyCalculator calc(1.0, 0.0);

  const Atom atom1(0,Vector3d(0,0,0), 0);
  std::vector<Atom> atoms1(1, atom1);
  const Molecule mol1(atoms1, "mol1", "testtest");

  const Atom atom2(0,Vector3d(4,0,0), 0);
  std::vector<Atom> atoms2(1, atom2);
  const Molecule mol2(atoms2, "mol2", "testtest");
  fltype expected = -0.035579*std::exp(-0.16)
    + (-0.005156) * std::exp(-1.96)
    + (-0.035069) * 1;

  fltype actual = calc.getEnergy(mol1, mol2);
  BOOST_CHECK_CLOSE(expected, actual, 2e-5);
}
