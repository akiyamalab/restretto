#include <istream>
#include <vector>
#include <memory>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/graphsym.h>
#include <openbabel/query.h>

#include "common.hpp"

#ifndef RMSD_H_
#define RMSD_H_

namespace OpenBabel {
	class AtomDistanceSorter {
    vector3 ref;
  public:
    AtomDistanceSorter(OpenBabel::OBAtom *r) : ref(r->GetVector()) {}
    bool operator()(OpenBabel::OBAtom *l, OpenBabel::OBAtom *r) const;
  };
  /* class to perform graph matching between two molecules.
  * Is initialized with the reference molecule.
  * Will figure out the atom correspondences and compute the rmsd between the
  * ref mol and a provided test mol.
  */
  class Matcher {
    const OpenBabel::OBMol *ref;
    std::shared_ptr<OpenBabel::OBQuery> query;
    std::shared_ptr<OpenBabel::OBIsomorphismMapper> mapper;

    class MapRMSDFunctor: public OpenBabel::OBIsomorphismMapper::Functor {
    private:
      const OpenBabel::OBMol *ref;
      OpenBabel::OBMol& test;
      double bestRMSD;
      bool minimize;
    public:
      //will modify t if min is true
      MapRMSDFunctor(const OpenBabel::OBMol* r, OpenBabel::OBMol& t, bool min = false) : ref(r), test(t), bestRMSD(HUGE_VAL), minimize(min) {}

      bool operator()(OpenBabel::OBIsomorphismMapper::Mapping &map) override;

      double getRMSD() const { return bestRMSD; }
    };

  public:
    Matcher(OpenBabel::OBMol& mol) : ref(&mol) {
      query = std::shared_ptr<OpenBabel::OBQuery>(CompileMoleculeQuery(&mol));
      mapper = std::shared_ptr<OpenBabel::OBIsomorphismMapper>(OpenBabel::OBIsomorphismMapper::GetInstance(query.get()));
    }

    //computes a correspondence between the ref mol and test (exhaustively)
    //and returns the rmsd; returns infinity if unmatchable
    double computeRMSD(OpenBabel::OBMol& test, bool minimize = false) {
      MapRMSDFunctor funct(ref, test, minimize);
      mapper->MapGeneric(funct, &test);
      return funct.getRMSD();
    }
  };
	
  // calculate minimum RMSD between mol and ref_mols
  fltype calc_minRMSD(const OpenBabel::OBMol& mol, const std::vector<OpenBabel::OBMol>& ref_mols);
}

#endif
