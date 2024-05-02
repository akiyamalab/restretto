#include <iostream>
#include <vector>
#include <iterator>
#include <string>
#include <sstream> 
#include <algorithm> //std::find in decomposite

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>

#include "Fragment.hpp"
#include "MoleculeToFragments.hpp"
#include "log_writer_stream.hpp"

#include "infile_reader.hpp"
#include "utils.hpp"

#include "OBMol.hpp"
#include "MoleculeToFragments.hpp"

namespace {
  struct DecomposeConfiguration {
    std::vector<std::string> ligand_files;
    std::string fragment_file, output_file, log_file;
    int capping_atomic_num, max_ring_size;
    bool do_carbon_capping, insert_fragment_id_to_isotope, merge_solitary;
  };

  DecomposeConfiguration parseArgs(int argc, char **argv){

    // Definition of options
    using namespace boost::program_options;
    options_description options("Options");
    options_description hidden("Hidden options");
    positional_options_description pos_desc;
    options.add_options()
      ("help,h", "show help")
      ("fragment,f", value<std::string>(), "fragment file (.mol2 file)")
      ("ligand,l", value<std::vector<std::string> >()->multitoken(), "ligand file")
      ("output,o", value<std::string>(), "output (annotated) ligand file (.sdf file)")
      ("log", value<std::string>(), "log file")
      ("capping_atomic_num", value<int>()->default_value(-1), "The atomic number of capping atoms. No capping if it is set to -1")
      ("enable_carbon_capping", bool_switch()->default_value(false), "Enabling capping even for Carbon atoms")
      ("ins_fragment_id", bool_switch()->default_value(false), "Enabling isotope number injection to mark fragment IDs")
      ("max_ring_size", value<int>()->default_value(-1), "Maximum ring size")
      ("no_merge_solitary", bool_switch()->default_value(false), "Disabling merging of solitary fragments");

    options_description desc;
    desc.add(options).add(hidden);
    variables_map vmap;
    store(command_line_parser(argc, argv).
          options(desc).positional(pos_desc).run(), vmap);
    notify(vmap);


    // showing help dialog
    bool insufficient_input = (!vmap.count("fragment") || !vmap.count("ligand") || !vmap.count("output"));
    bool help_mode = vmap.count("help");
    if (insufficient_input || help_mode){
      if (!help_mode){
        std::cout << "too few arguments" << std::endl;
      }
      
      std::cout << "Usage: " << argv[0] << "[options]" << std::endl
                << options << std::endl;
      std::exit(0);
    }
    
    // parse input options and configuration file
    DecomposeConfiguration conf;
    if (vmap.count("ligand")) conf.ligand_files = vmap["ligand"].as<std::vector<std::string> >();
    if (vmap.count("fragment")) conf.fragment_file = vmap["fragment"].as<std::string>();
    if (vmap.count("output")) conf.output_file = vmap["output"].as<std::string>();
    if (vmap.count("log")) conf.log_file = vmap["log"].as<std::string>();
    conf.capping_atomic_num = vmap["capping_atomic_num"].as<int>();
    conf.do_carbon_capping  = vmap["enable_carbon_capping"].as<bool>();
    conf.insert_fragment_id_to_isotope = vmap["ins_fragment_id"].as<bool>();
    conf.max_ring_size = vmap["max_ring_size"].as<int>();
    conf.merge_solitary = !vmap["no_merge_solitary"].as<bool>();

    return conf;
  }

  fragdock::Molecule convert_molecule(OpenBabel::OBMol& obmol) {
    obmol.AddHydrogens();
    fragdock::Molecule ligand_mol = format::toFragmentMol(obmol);
    obmol.DeleteHydrogens();
    // ligand_mol.translate(-ligand_mol.getCenter());
    return ligand_mol;
  }

  void correctvalence(OpenBabel::OBMol* obmol){
    OpenBabel::OBConversion conv_mol2can, conv_can2mol;
    conv_mol2can.SetInFormat("mol");
    conv_can2mol.SetInFormat("can");
    conv_mol2can.SetOutFormat("can");
    conv_can2mol.SetOutFormat("mol");
    std::string mol_string = conv_can2mol.WriteString(obmol);
    conv_mol2can.ReadString(obmol, mol_string);
  }

  int decomposite(const std::vector<OpenBabel::OBMol>&  molecules,
                        std::vector<std::string>&       frag_smi_list,
                        std::vector<OpenBabel::OBMol>&  annotated_mols,
                        std::vector<OpenBabel::OBMol>&  fragments,
                        int                             capping_atomic_num = -1,
                        bool                            capping_for_carbon = false,
                        int                             max_ring_size = -1,
                        bool                            merge_solitary = true,
                        bool                            insert_fragment_id_to_isotope=false){

    // TODO Segmentation fault was occurred when this function was parallelized.
    for(int i=0; i<molecules.size(); ++i){
      std::vector<std::string> frag_smiles;
      OpenBabel::OBMol obmol = molecules[i];
      fragdock::Molecule mol = convert_molecule(obmol);
      std::vector<fragdock::Fragment> frags = fragdock::DecomposeMolecule(mol, max_ring_size, merge_solitary, true);
      
      for(std::vector<fragdock::Fragment>::iterator fit=frags.begin(); fit!=frags.end(); ++fit){
        OpenBabel::OBMol obmol = format::toOBMol(*fit, molecules[i], capping_atomic_num, capping_for_carbon, insert_fragment_id_to_isotope);
        correctvalence(&obmol);
        std::string smiles = OpenBabel::canonicalSmiles(obmol);
        frag_smiles.push_back(smiles);
        obmol.SetTitle(smiles);
        
        if(std::find(frag_smi_list.begin(), frag_smi_list.end(), smiles) == frag_smi_list.end()){
          frag_smi_list.push_back(smiles);
          fragments.push_back(obmol);
        }

        int frag_num = (std::find(frag_smi_list.begin(), frag_smi_list.end(), smiles) - frag_smi_list.begin());
      }

      // add fragment info to molecules
      obmol = molecules[i];
      OpenBabel::SetProperty(obmol, "fragment_info", boost::algorithm::join(frag_smiles, ","));
      annotated_mols.push_back(obmol);

      // show progress
      if((i+1)%1000==0){
        logs::lout << logs::info << (i+1) << " compounds were finished to decompose. The number of fragments is " << fragments.size() << std::endl;
      }
    }

    if(molecules.size()%1000!=0){
      logs::lout << logs::info << molecules.size() << " compounds were finished to decompose. The number of fragments is " << fragments.size() << std::endl;
    }

    return fragments.size();
  }

  void outputmolecules(const std::vector<OpenBabel::OBMol>& molecules,
                       const std::string& filename){
    OpenBabel::outputOBMol outputs(filename);
    for (int i=0; i<molecules.size(); i++) {
      outputs.write(molecules[i]);
    }
    outputs.close();
  }

}


int main(int argc, char** argv){

  // config load
  DecomposeConfiguration config = parseArgs(argc, argv);

  // logging start
  if(config.log_file == ""){
    config.log_file = (config.output_file + "__" + utils::getDate() + ".log");
  }
  logs::log_init(config.log_file, true);

  // preparing ligand data
  std::vector<std::vector<int> > mol2frag_list;
  std::vector<OpenBabel::OBMol> fragments;
  std::vector<std::string> frag_smi_list;
  std::vector<OpenBabel::OBMol> annotated_mols;
  for (int i=0; i<config.ligand_files.size(); i++){
    logs::lout << logs::info << "parse ligand file: " << config.ligand_files[i] << std::endl;
  }
  std::vector<OpenBabel::OBMol> molecules = format::ParseFileToOBMol(config.ligand_files);
  logs::lout << logs::info << "decomposite ligands into fragments" << std::endl;
  decomposite(molecules, frag_smi_list, annotated_mols, fragments, 
              config.capping_atomic_num, config.do_carbon_capping, config.max_ring_size, 
              config.merge_solitary, config.insert_fragment_id_to_isotope); // TODO: use merge_solitary

  logs::lout << logs::info << "output ligands added fragment information at " << config.output_file << std::endl;
  outputmolecules(annotated_mols, config.output_file);

  logs::lout << logs::info << "output fragments at " << config.fragment_file << std::endl;
  outputmolecules(fragments,      config.fragment_file);

  logs::close;
  return 0;
}