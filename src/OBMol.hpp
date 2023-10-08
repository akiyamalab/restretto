#include <istream>
#include <vector>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/graphsym.h>
#include <openbabel/query.h>

#include "Molecule.hpp"

#ifndef OBMOL_H_
#define OBMOL_H_

namespace format{
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(const std::vector<std::string>& filepaths);
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(const std::string& filepath);
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(std::istream& stream, const std::string& in_format);
  fragdock::Molecule    toFragmentMol(const OpenBabel::OBMol& mol);
  OpenBabel::OBMol      toOBMol(const fragdock::Molecule &mol, const OpenBabel::OBMol& original_obmol);
}

namespace OpenBabel{
  class outputOBMol {
    std::ofstream ofs;
    OpenBabel::OBConversion conv;
  public:
    outputOBMol(const std::string& filepath) {
      ofs.open(filepath);
      conv = OpenBabel::OBConversion(NULL, &ofs);
      conv.SetOutFormat("sdf");
    }
    void write(const OpenBabel::OBMol& mol) {
      conv.Write(&const_cast<OpenBabel::OBMol&>(mol));
    }
    void close() {
      ofs.close();
    }
  };
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
  std::string canonicalSmiles(const OpenBabel::OBMol& mol);
  void toCanonical(OpenBabel::OBMol& mol);
  std::vector<uint> getRenumber(OpenBabel::OBMol& mol);
  void UpdateCoords(OpenBabel::OBMol& mol, const fragdock::Molecule& ref_mol);
  void outputOBMolsToSDF(const std::string& filepath, const std::vector<OpenBabel::OBMol>& mols);
  void outputOBMolsToSDF(const std::string& filepath, const std::vector<std::vector<OpenBabel::OBMol> >& molss);
  OpenBabel::OBMol CreateSubMolecule(const OpenBabel::OBMol& reference, const std::vector<int>& atom_ids);
  //preprocess molecule into a standardized state for heavy atom rmsd computation
  void processMol(OpenBabel::OBMol& mol);
}

#endif
