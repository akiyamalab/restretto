#include "common.hpp"
#include "utils.hpp"
#include "OBMol.hpp"
#include "infile_reader.hpp"
#include "EnergyCalculator.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <sys/stat.h>
#include <chrono>

namespace {
  format::DockingConfiguration parseArgs(int argc, char **argv){
    using namespace boost::program_options;
    options_description options("Options");
    options_description hidden("Hidden options");
    positional_options_description pos_desc;
    hidden.add_options()
      ("conf-file", value<std::string>(), "configuration file");
    pos_desc.add("conf-file", 1);
    options.add_options()
      ("help,h", "show help")
      ("ligand,l", value<std::vector<std::string> >()->multitoken(), "ligand file");
      // ("receptor,r", value<std::string>(), "receptor file (.pdb file)");
    options_description desc;
    desc.add(options).add(hidden);
    variables_map vmap;
    store(command_line_parser(argc, argv).
	  options(desc).positional(pos_desc).run(), vmap);
    notify(vmap);

    if (!vmap.count("conf-file") || vmap.count("help")){
      if (!vmap.count("conf-file") && !vmap.count("help")){
	std::cout << "too few arguments" << std::endl;
      }
      std::cout << "Usage: ligandock conf-file [options]\n"
		<< options << std::endl;
      std::exit(0);
    }
    format::DockingConfiguration conf = format::ParseInFile(vmap["conf-file"].as<std::string>().c_str());
    if (vmap.count("ligand")) conf.ligand_files = vmap["ligand"].as<std::vector<std::string> >();
    // if (vmap.count("receptor")) conf.receptor_file = vmap["receptor"].as<std::string>();
    return conf;
  }

} // namespace

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;

  format::DockingConfiguration config = parseArgs(argc, argv);


  // // parse receptor file
  // OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  // Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse ligands file
  vector<OpenBabel::OBMol> ligands = format::ParseFileToOBMol(config.ligand_files);

  int ligs_sz = ligands.size();


  EnergyCalculator calc;
  // EnergyCalculator calc(1.0, 0.0);

  cout << fixed << setprecision(5);

  for (int lig_ind = 0; lig_ind < ligs_sz; ++lig_ind) {
    OpenBabel::OBMol& ligand = ligands[lig_ind];
    ligand.AddHydrogens();
    Molecule mol = format::toFragmentMol(ligand);
    ligand.DeleteHydrogens();
    mol.deleteHydrogens();

    cout << "Title: " << mol.gettitle() << '\n';
    cout << "IntraEnergy: " << EnergyCalculator::getIntraEnergy_strict(mol) << '\n';
  }

  return 0;
}
