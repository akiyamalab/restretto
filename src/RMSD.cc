#include "RMSD.hpp"

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
        //with minimize on, change coordinates
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

  void processMol(OpenBabel::OBMol& mol) {
    //isomorphismmapper wants isomorphic atoms to have the same aromatic and ring state,
    //but these proporties aren't reliable enough to be trusted in evaluating molecules
    //should be considered the same based solely on connectivity
    
    mol.DeleteHydrogens(); //heavy atom rmsd
    for(OpenBabel::OBAtomIterator aitr = mol.BeginAtoms(); aitr != mol.EndAtoms(); aitr++) {
      OpenBabel::OBAtom *a = *aitr;
      a->UnsetAromatic();
      a->SetInRing();
    }
    for(OpenBabel::OBBondIterator bitr = mol.BeginBonds(); bitr != mol.EndBonds(); bitr++) {
      OpenBabel::OBBond *b = *bitr;
      b->UnsetAromatic();
      b->SetBondOrder(1);
      b->SetInRing();
    }
    //avoid recomputations
    mol.SetHybridizationPerceived();
    mol.SetRingAtomsAndBondsPerceived();
    mol.SetAromaticPerceived();
  }
}