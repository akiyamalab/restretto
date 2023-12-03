#include <istream>

#include <openbabel/canon.h>

#include "OBMol.hpp"
#include "utils.hpp"

namespace {
  int getXSType(const OpenBabel::OBAtom& oatom) {
    OpenBabel::OBAtom& atom = const_cast<OpenBabel::OBAtom&>(oatom);

    if      (atom.IsCarbon()   && !atom.GetHeteroValence()) return XS_TYPE_C_H;
    else if (atom.IsCarbon()   &&  atom.GetHeteroValence()) return XS_TYPE_C_P;
    else if (atom.IsNitrogen() && !atom.IsHbondAcceptor() && !atom.IsHbondDonor()) return XS_TYPE_N_P;
    else if (atom.IsNitrogen() && !atom.IsHbondAcceptor() &&  atom.IsHbondDonor()) return XS_TYPE_N_D;
    else if (atom.IsNitrogen() &&  atom.IsHbondAcceptor() && !atom.IsHbondDonor()) return XS_TYPE_N_A;
    else if (atom.IsNitrogen() &&  atom.IsHbondAcceptor() &&  atom.IsHbondDonor()) return XS_TYPE_N_DA;
    else if (atom.IsOxygen()   && !atom.IsHbondAcceptor() && !atom.IsHbondDonor()) return XS_TYPE_O_P;
    else if (atom.IsOxygen()   && !atom.IsHbondAcceptor() &&  atom.IsHbondDonor()) return XS_TYPE_O_D;
    else if (atom.IsOxygen()   &&  atom.IsHbondAcceptor() && !atom.IsHbondDonor()) return XS_TYPE_O_A;
    else if (atom.IsOxygen()   &&  atom.IsHbondAcceptor() &&  atom.IsHbondDonor()) return XS_TYPE_O_DA;
    else if (atom.IsSulfur())           return XS_TYPE_S_P;
    else if (atom.IsPhosphorus())       return XS_TYPE_P_P;
    else if (atom.GetAtomicNum() ==  9) return XS_TYPE_F_H;
    else if (atom.GetAtomicNum() == 17) return XS_TYPE_Cl_H;
    else if (atom.GetAtomicNum() == 35) return XS_TYPE_Br_H;
    else if (atom.GetAtomicNum() == 53) return XS_TYPE_I_H;
    else if (atom.IsMetal())            return XS_TYPE_Met_D;
    else if (atom.IsHydrogen())         return XS_TYPE_H;
    else {
      logs::lout << "Unknown atom is read." << std::endl;
      logs::lout << "Atomic Number: " << atom.GetAtomicNum() << std::endl;
      logs::lout << "Atomic Type: " << atom.GetType() << std::endl;
      return XS_TYPE_OTHER;
    }
    // else return XS_TYPE_SIZE;
  }

}

namespace format{
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(const std::vector<std::string>& filepaths) {
    std::vector<OpenBabel::OBMol> ret;
    for (auto& path : filepaths) {
      std::vector<OpenBabel::OBMol> mols = ParseFileToOBMol(path);
      ret.insert(ret.end(), mols.begin(), mols.end());
    }
    return ret;
  }
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(const std::string& filepath){
    std::ifstream ifs(filepath.c_str());
    if(ifs.fail()){
      std::cerr << "failed to open: " << filepath << std::endl;
      abort();
    }
    std::string format = utils::GetExtension(filepath);
    return ParseFileToOBMol(ifs, format);
  }
  std::vector<OpenBabel::OBMol> ParseFileToOBMol(std::istream& stream, const std::string& in_format){
    OpenBabel::OBConversion conv(&stream);
    conv.SetInFormat(in_format.c_str());
    OpenBabel::OBMol mol;
    std::vector<OpenBabel::OBMol> molecules;
    for(bool end_frag = !conv.Read(&mol); end_frag != true; end_frag = !conv.Read(&mol)){
      mol.AddPolarHydrogens();
      molecules.push_back(mol);
    }

    return molecules;
  }

  fragdock::Molecule toFragmentMol(const OpenBabel::OBMol &c_obmol){
    OpenBabel::OBMol &obmol = const_cast<OpenBabel::OBMol&>(c_obmol);
    std::vector<fragdock::Atom> fatoms(obmol.NumAtoms());


    for(int i=0; i<obmol.NumAtoms(); i++) {
      OpenBabel::OBAtom& oatom = *( obmol.GetAtom(i+1) ); //obmol::index is 1 origin
      int xstype = getXSType(oatom);

      fragdock::Vector3d pos(oatom.GetX(), oatom.GetY(), oatom.GetZ());
      assert(i == oatom.GetId());
      fatoms[i] = fragdock::Atom(oatom.GetId(), pos, xstype);
    }

    fragdock::Molecule mol(fatoms, obmol.GetTitle(), OpenBabel::canonicalSmiles(obmol));
    for(int i = 0; i < obmol.NumBonds(); i++){
      mol.append(fragdock::Bond(obmol.GetBond(i)->GetBeginAtom()->GetId(),
				obmol.GetBond(i)->GetEndAtom()->GetId(),
				obmol.GetBond(i)->IsRotor() ));
				// obmol.GetBond(i)->IsRotor() || (obmol.GetBond(i)->IsSingle() && !obmol.GetBond(i)->IsInRing()) ));
    }

    return mol;
  }

  OpenBabel::OBMol toOBMol(const fragdock::Molecule &mol, const OpenBabel::OBMol& original_obmol) {
    using namespace fragdock;
    using namespace std;

    OpenBabel::OBMol obmol = original_obmol;
    obmol.SetTitle(mol.gettitle().c_str());

    std::vector<int> del_id_list(original_obmol.NumAtoms());
    for(int i=0; i<original_obmol.NumAtoms(); i++){
      del_id_list[i] = i;
    }


    for(int i=0; i<mol.getAtoms().size(); i++) {
      int id = mol.getAtoms()[i].getId();
      if (mol.getAtoms()[i].getXSType() != XS_TYPE_H)
        del_id_list.erase(std::remove(del_id_list.begin(), del_id_list.end(), id), del_id_list.end());

      if (mol.getAtoms()[i].getXSType() == XS_TYPE_DUMMY) {
        obmol.GetAtom(id + 1)->SetFormalCharge(0);
        obmol.GetAtom(id + 1)->SetAtomicNum(0);
      }
      else if (mol.getAtoms()[i].getXSType() != XS_TYPE_H) {
        obmol.GetAtom(id + 1)->SetIsotope(mol.getAtoms()[i].getXSType());
      }
    }

    for(int i=0; i<del_id_list.size(); i++){
      int id = original_obmol.GetAtom(del_id_list[i] + 1)->GetId();
      obmol.DeleteAtom(obmol.GetAtomById(id));
    }

    /*
    vector<Atom> atoms = mol.getAtoms();
    for(int i = 0; i < obmol.NumAtoms(); i++) {
      OpenBabel::OBAtom* atom = obmol.GetAtomById(i);
      double x = atoms[i].x;
      double y = atoms[i].y;
      double z = atoms[i].z;
      atom->SetVector(x, y, z);
    }
    */


    return obmol;
  }

}

namespace OpenBabel{
  std::string canonicalSmiles(const OpenBabel::OBMol& mol){
    OpenBabel::OBConversion conv;
    conv.SetOutFormat("can");
    conv.AddOption("n");
    return utils::trim_tail(conv.WriteString(&const_cast<OpenBabel::OBMol&>(mol)), " \t\n");
  }

  void toCanonical(OpenBabel::OBMol& mol) {
    OpenBabel::OBMol * pmol = &mol;
    std::vector<OpenBabel::OBAtom*> atoms;
    FOR_ATOMS_OF_MOL (atom, pmol)
      atoms.push_back(&*atom);

    std::vector<unsigned int> symmetry_classes;
    OpenBabel::OBGraphSym gs(pmol);
    gs.GetSymmetry(symmetry_classes);

    std::vector<unsigned int> canon_labels;
    CanonicalLabels(pmol, symmetry_classes, canon_labels);

    std::vector<OpenBabel::OBAtom*> newatoms(atoms.size(), 0);
    for (int i = 0; i < canon_labels.size(); ++i)
      newatoms[canon_labels[i] - 1] = atoms[i];

    pmol->RenumberAtoms(newatoms);
  }

  /**
   * Determine the canonical order indices of atoms in the molecule.
   * @param[in] mol the molecule
   * @return The canonical order indices of atoms in the molecule
  */
  std::vector<uint> getRenumber(OpenBabel::OBMol& mol) {
    std::vector<uint> canon_labels;
    OpenBabel::OBMol * pmol = &mol;
    std::vector<unsigned int> symmetry_classes;
    OpenBabel::OBGraphSym gs(pmol);
    gs.GetSymmetry(symmetry_classes);

    OpenBabel::CanonicalLabels(pmol, symmetry_classes, canon_labels);
    return canon_labels;
  }

  void UpdateCoords(OpenBabel::OBMol& mol, const fragdock::Molecule& ref_mol){
    const std::vector<fragdock::Atom>& ref_atoms = ref_mol.getAtoms();
    for(int i=0; i<mol.NumAtoms(); i++){
      const fragdock::Atom& ref_atom = ref_atoms[i];
      mol.GetAtom(i + 1)->SetVector(ref_atom.x, ref_atom.y, ref_atom.z);
    }
  }

  void outputOBMolsToSDF(const std::string& filepath, const std::vector<OpenBabel::OBMol>& mols) {
    std::ofstream output;
    output.open(filepath);
    OpenBabel::OBConversion conv(NULL, &output);
    conv.SetOutFormat("sdf");
    for(const OpenBabel::OBMol& mol: mols){
      conv.Write(&const_cast<OpenBabel::OBMol&>(mol));
    }
    output.close();
  }
  void outputOBMolsToSDF(const std::string& filepath, const std::vector<std::vector<OpenBabel::OBMol> >& molss) {
    std::ofstream output;
    output.open(filepath);
    OpenBabel::OBConversion conv(NULL, &output);
    conv.SetOutFormat("sdf");
    for (auto& mols : molss) {
      for (auto& mol : mols) {
        conv.Write(&const_cast<OpenBabel::OBMol&>(mol));
      }
    }
    output.close();
  }

  OpenBabel::OBMol CreateSubMolecule(const OpenBabel::OBMol& reference, const std::vector<int>& atom_ids){
    using namespace std;

    OpenBabel::OBMol mol = reference;

    std::set<int> del_ids;
    for(int i=0; i<reference.NumAtoms(); i++) del_ids.insert(i);

    for(int elem: atom_ids){
      del_ids.erase(elem);
    }

    for(int elem: del_ids){
      int id = reference.GetAtom(elem + 1)->GetId();
      mol.DeleteAtom(mol.GetAtomById(id));
    }
    return mol;
  }
  void SetProperty(OpenBabel::OBMol& mol, const std::string& key, const std::string& value){
    OpenBabel::OBPairData *dp = new OpenBabel::OBPairData;
    dp->SetAttribute(key);
    dp->SetValue(value);
    mol.SetData(dp);
  }
  void SetProperty(OpenBabel::OBMol& mol, const std::string& key, fltype value){
    SetProperty(mol, key, std::to_string(value));
  }
}
