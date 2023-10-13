#include "OBMol.hpp"
#include "common.hpp"
#include "infile_reader.hpp"
#include "utils.hpp"

#include "AtomInterEnergyGrid.hpp"
#include "EnergyCalculator.hpp"
#include "log_writer_stream.hpp"

#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <sys/stat.h>

#include <boost/program_options.hpp>
#include <boost/format.hpp>

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
      ("grid,g", value<std::string>(), "grid folder")
      ("receptor,r", value<std::string>(), "receptor file (.pdb file)");
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
      std::cout << "Usage: ligangrid conf-file -o folderName [options]\n"
		<< options << std::endl;
      std::exit(0);
    }

    format::DockingConfiguration conf = format::ParseInFile(vmap["conf-file"].as<std::string>().c_str());
    if (vmap.count("grid")) conf.grid_folder = vmap["grid"].as<std::string>();
    if (vmap.count("receptor")) conf.receptor_file = vmap["receptor"].as<std::string>();
    return conf;
  }

  void makeFolder(std::string folderName){
    struct stat st;
    if(stat(folderName.c_str(), &st)==-1){
      std::cout << "there isn't folder. Try create: " << folderName << std::endl;
      if(mkdir(folderName.c_str(), 0775) == -1){
	std::cerr << "Fail to make folder: " << folderName << std::endl;
	abort();
      }
    }else if(!S_ISDIR(st.st_mode)){
      std::cerr << "Fail to make folder. There is same name file (not folder): " << folderName << std::endl;
      abort();
    }
  }

  void makeEnergyGrid(fragdock::AtomInterEnergyGrid& energy_grid, const fragdock::Molecule& receptor_mol, const fragdock::EnergyCalculator& calc){
    const fragdock::Vector3d temp(0,0,0);
    fragdock::Atom atom(0, temp, energy_grid.getXSType());

    // #pragma omp parallel for // TODO: the result is changed when parallelized
    for(int x=0; x<energy_grid.getNum().x; x++) {
      for(int y=0; y<energy_grid.getNum().y; y++) {
        for(int z=0; z<energy_grid.getNum().z; z++) {
          atom.setPos(energy_grid.convert(x, y, z));
          energy_grid.setEnergy(x, y, z, calc.getEnergy(atom, receptor_mol));
        }
      }
    }

  }

} // namespace


int main(int argc, char **argv){
  const format::DockingConfiguration conf = parseArgs(argc, argv);
  makeFolder(conf.grid_folder);

  //Prepareing logger
  logs::log_init("atomgrid-gen.log");

  // logs::lout << logs::info << "Read parameter file. " << conf.forcefield_file;
  // fragdock::AtomType::setAtomParams(conf.forcefield_file.c_str());

  logs::lout << logs::info << "Prepare receptor molecule.";
  const OpenBabel::OBMol receptor = format::ParseFileToOBMol(conf.receptor_file.c_str())[0];
  const fragdock::Molecule receptor_mol = format::toFragmentMol(receptor);

  const format::SearchGrid& grid = conf.grid;

  const fragdock::EnergyCalculator calc(0.95);

  //Making energyGrid
  #pragma omp parallel for // making energy grids can be parallelized
  for(int xs_type = 0; xs_type < XS_TYPE_SIZE; xs_type++){
    logs::lout << logs::info << "Prepare energy grid. " << xs_strings[xs_type] << std::endl;

    //prepare energy_grid
    // const fltype MARGIN = 10;
    // const fltype DEFAULT_PITCH = 0.25; // (atomtypes[i]=="e") ? 0.5 : 0.25;
    const fragdock::Point3d<fltype>& center = grid.center;
    const fragdock::Point3d<fltype>& pitch = grid.score_pitch;
    const fragdock::Point3d<int> num = utils::ceili(grid.outer_width / 2 / pitch) * 2 + 1;

    logs::lout << logs::debug << "Initialize fragdock::AtomInterEnergyGrid. " << xs_strings[xs_type] << std::endl;
    fragdock::AtomInterEnergyGrid energy_grid(center, pitch, num, xs_type);
    logs::lout << logs::debug << "Start calculation of each point. " << xs_strings[xs_type] << std::endl;
    makeEnergyGrid(energy_grid, receptor_mol, calc);

    logs::lout << logs::debug << "Make energy grid file. " << xs_strings[xs_type] << std::endl;
    energy_grid.writeFile(conf.grid_folder + "/" + xs_strings[xs_type] + ".grid");
    logs::lout << logs::info << "Making grid file has done. " << xs_strings[xs_type] << std::endl;
  }

  return 0;
}
