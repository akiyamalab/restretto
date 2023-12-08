#include "RMSD.hpp"

/*
 Reference: openbabel-2-4, openbabel/tools/orbrms.cpp
 GitHub URL: https://github.com/openbabel/openbabel/blob/openbabel-2-4-x/tools/obrms.cpp
*/

namespace OpenBabel{
  bool AtomDistanceSorter::operator()(OpenBabel::OBAtom *l, OpenBabel::OBAtom *r) const {
    double ld = ref.distSq(l->GetVector());
    double rd = ref.distSq(r->GetVector());

    return ld < rd;
  }

  bool OpenBabel::Matcher::MapRMSDFunctor::operator()(OpenBabel::OBIsomorphismMapper::Mapping &map) {
    unsigned N = map.size();
    double *refcoord = (double*)alloca(sizeof(double)*N * 3);
    double *testcoord = (double*)alloca(sizeof(double)*N * 3);

    for (unsigned i = 0; i < N; i++) {
      //obmol indices are 1-indexed while the mapper is zero indexed
      const OBAtom *ratom = ref->GetAtom(map[i].first + 1);
      const OBAtom *tatom = test.GetAtom(map[i].second + 1);
      assert(ratom && tatom);

      for (unsigned c = 0; c < 3; c++) {
        refcoord[3 * i + c] = ratom->GetVector()[c];
        testcoord[3 * i + c] = tatom->GetVector()[c];
      }
    }

    if (minimize) {
      double rmatrix[3][3] = { 0 };

      double rave[3] = { 0, 0, 0 };
      double tave[3] = { 0, 0, 0 };
      //center
      for (unsigned i = 0; i < N; i++) {
        for (unsigned c = 0; c < 3; c++) {
          rave[c] += refcoord[3 * i + c];
          tave[c] += testcoord[3 * i + c];
        }
      }

      for (unsigned c = 0; c < 3; c++) {
        rave[c] /= N;
        tave[c] /= N;
      }

      for (unsigned i = 0; i < N; i++) {
        for (unsigned c = 0; c < 3; c++) {
          refcoord[3 * i + c] -= rave[c];
          testcoord[3 * i + c] -= tave[c];
        }
      }

      qtrfit(refcoord, testcoord, N, rmatrix);
      rotate_coords(testcoord, rmatrix, N); 

      for (unsigned i = 0; i < N; i++) {
        //with minimize on, change coordinatescomputeRMSD
        OBAtom *tatom = test.GetAtom(map[i].second + 1);
        tatom->SetVector(testcoord[3*i]+rave[0], testcoord[3*i+1]+rave[1], testcoord[3*i+2]+rave[2]);
      }
    }

    double rmsd = calc_rms(refcoord, testcoord, N);

    if (rmsd < bestRMSD) {
      bestRMSD = rmsd;
    }
    // check all possible mappings
    return false;
  }

  OpenBabel::OBMol standardize_mol(const OpenBabel::OBMol& mol) {
    OpenBabel::OBMol m = mol; // copying object to avoid modifying the original

    /* isomorphismmapper wants isomorphic atoms to have the same aromatic and ring state,
     * but these proporties aren't reliable enough to be trusted in evaluating molecules
     * should be considered the same based solely on connectivity
    */
    m.DeleteHydrogens(); //heavy atom rmsd
    for(OpenBabel::OBAtomIterator aitr = m.BeginAtoms(); aitr != m.EndAtoms(); aitr++) {
      OpenBabel::OBAtom *a = *aitr;
      a->UnsetAromatic();
      a->SetInRing();
    }
    for(OpenBabel::OBBondIterator bitr = m.BeginBonds(); bitr != m.EndBonds(); bitr++) {
      OpenBabel::OBBond *b = *bitr;
      b->UnsetAromatic();
      b->SetBondOrder(1);
      b->SetInRing();
    }
    //avoid recomputations
    m.SetHybridizationPerceived();
    m.SetRingAtomsAndBondsPerceived();
    m.SetAromaticPerceived();

    return m;
  }

  fltype calc_minRMSD(const OpenBabel::OBMol& mol, const std::vector<OpenBabel::OBMol>& ref_mols) {

    OpenBabel::OBMol m = standardize_mol(mol);
    std::vector<OpenBabel::OBMol> cf_mols = ref_mols;
    for (OpenBabel::OBMol& cf_mol : cf_mols) {
      cf_mol = standardize_mol(cf_mol);
    }

    // calculate minimum RMSD between mol and ref_mols
    OpenBabel::Matcher matcher(m);
    fltype minRMSD = HUGE_VAL;
    for (OpenBabel::OBMol cf_mol : cf_mols) { 
      fltype rmsd = matcher.computeRMSD(cf_mol);
      logs::lout << "RMSD: " << rmsd << std::endl;
      minRMSD = std::min(minRMSD, rmsd);
    }
    return minRMSD;
  }
}