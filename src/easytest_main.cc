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

#include <unordered_map>

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
      ("output,o", value<std::string>(), "output file (.mol2 file)")
      ("ligand,l", value<std::vector<std::string> >()->multitoken(), "ligand file")
      ("receptor,r", value<std::string>(), "receptor file (.pdb file)")
      ("grid,g", value<std::string>(), "grid folder")
      ("memsize,m", value<int64_t>(), "fragment grid's memory size[MB]");
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
      std::cout << "Usage: ligandock conf-file [options]\n" << options << std::endl;
      std::exit(0);
    }
    format::DockingConfiguration conf = format::ParseInFile(vmap["conf-file"].as<std::string>().c_str());
    if (vmap.count("ligand")) conf.ligand_files = vmap["ligand"].as<std::vector<std::string> >();
    if (vmap.count("receptor")) conf.receptor_file = vmap["receptor"].as<std::string>();
    if (vmap.count("output")) conf.output_file = vmap["output"].as<std::string>();
    if (vmap.count("grid")) conf.grid_folder = vmap["grid"].as<std::string>();
    if (vmap.count("memsize")) conf.mem_size = vmap["memsize"].as<int64_t>();
    return conf;
  }

} // namespace

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;

  format::DockingConfiguration config = parseArgs(argc, argv);


  // parse receptor file
  OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse ligands file
  vector<OpenBabel::OBMol> ligands = format::ParseFileToOBMol(config.ligand_files);

  int ligs_sz = ligands.size();


  vector<Molecule> ligands_mol(ligs_sz);

  EnergyCalculator calc(1.0, 0.0);

  for (int lig_ind = 0; lig_ind < ligs_sz; ++lig_ind) {
    OpenBabel::OBMol& ligand = ligands[lig_ind];
    ligand.AddHydrogens();
    ligands_mol[lig_ind] = format::toFragmentMol(ligand);
    ligand.DeleteHydrogens();

    Molecule& mol = ligands_mol[lig_ind];

    mol.deleteHydrogens();
    mol.setIntraEnergy(calc.calcIntraEnergy(mol));

    // mol.deleteHydrogens();
    // mol.calcRadius();
    fltype score = calc.getEnergy(mol, receptor_mol);
    fltype final_score = score / (1 + 0.05846 * mol.getNrots());
    cout << final_score << '\n';
    Vector3d center = mol.getCenter();
    // mol.translate(-ofst);
    const fltype pi = acos(-1.0);
    const int NEAREST_NUM = 200;
    fltype trans_step = 0.5;
    fltype rotate_step = pi / 30.0;
    for(int i = 0; i < NEAREST_NUM; i++) {
      Vector3d dv = Vector3d(utils::randf(-trans_step, trans_step),
                              utils::randf(-trans_step, trans_step),
                              utils::randf(-trans_step, trans_step));
      fltype th  = utils::randf(-rotate_step, rotate_step);
      fltype phi = utils::randf(-rotate_step, rotate_step);
      fltype psi = utils::randf(-rotate_step, rotate_step);

      Molecule tmp_mol = mol;
      tmp_mol.translate(-center);
      tmp_mol.rotate(th, phi, psi);
      tmp_mol.translate(center + dv);
      // tmp_mol.translate(dv);

      fltype val = calc.getEnergy(tmp_mol, receptor_mol);
      fltype final_val = val / (1 + 0.05846 * tmp_mol.getNrots());
      cout << final_val << '\n';
    }
  }

  return 0;
}
