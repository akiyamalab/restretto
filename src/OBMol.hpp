#include <istream>
#include <vector>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/babelconfig.h>
#include <openbabel/op.h>
#include <openbabel/graphsym.h>

#include "Molecule.hpp"
#include "utils.hpp"

#ifndef OBMOL_H_
#define OBMOL_H_

namespace format{
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(const std::vector<std::string>& filepaths);
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(const std::string& filepath);
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(std::istream& stream, const std::string& in_format);
  fragdock::Molecule    toFragmentMol(const OpenBabel::OBMol& mol);
  OpenBabel::OBMol      toOBMol(const fragdock::Molecule &mol, 
                                const OpenBabel::OBMol& original_obmol, 
                                int capping_atomic_num=-1, 
                                bool capping_for_carbon=false,
                                bool insert_fragment_id_to_isotope=true);
}

namespace OpenBabel{
  class outputOBMol {
    std::ofstream ofs;
    OpenBabel::OBConversion conv;
  public:
    outputOBMol(const std::string& filepath) {
      ofs.open(filepath);
      conv = OpenBabel::OBConversion(NULL, &ofs);
      // conv.SetOutFormat("sdf");
      conv.SetOutFormat(utils::GetExtension(filepath).c_str());
    }
    void write(const OpenBabel::OBMol& mol) {
      conv.Write(&const_cast<OpenBabel::OBMol&>(mol));
    }
    void close() {
      ofs.close();
    }
  };
  std::string canonicalSmiles(const OpenBabel::OBMol& mol);
  void toCanonical(OpenBabel::OBMol& mol);
  std::vector<uint> getRenumber(OpenBabel::OBMol& mol);
  void UpdateCoords(OpenBabel::OBMol& mol, const fragdock::Molecule& ref_mol);
  void outputOBMolsToSDF(const std::string& filepath, const std::vector<OpenBabel::OBMol>& mols);
  void outputOBMolsToSDF(const std::string& filepath, const std::vector<std::vector<OpenBabel::OBMol> >& molss);
  OpenBabel::OBMol CreateSubMolecule(const OpenBabel::OBMol& reference, const std::vector<int>& atom_ids);
  void SetProperty(OpenBabel::OBMol& mol, const std::string& key, const std::string& value);
  void SetProperty(OpenBabel::OBMol& mol, const std::string& key, fltype value);
}

#endif
