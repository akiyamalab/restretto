#include "common.hpp"
#include "utils.hpp"
#include "OBMol.hpp"
#include "infile_reader.hpp"
#include "log_writer_stream.hpp"
#include "AtomEnergyGrid.hpp"
#include "EnergyCalculator.hpp"
#include "Optimizer.hpp"

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
      ("grid,g", value<std::string>(), "grid folder");
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
    if (vmap.count("receptor")) conf.receptor_file = vmap["receptor"].as<std::string>();
    if (vmap.count("output")) conf.output_file = vmap["output"].as<std::string>();
    if (vmap.count("grid")) conf.grid_folder = vmap["grid"].as<std::string>();
    return conf;
  }

  std::string getDate() {
    time_t timer = time(NULL);
    tm* date = localtime(&timer);
    std::string dateStr = (boost::format("%02d_%02d_%02d_%02d_%02d")
			   % date->tm_mon % date->tm_mday
			   % date->tm_hour % date->tm_min % date-> tm_sec).str();
    return dateStr;
  }

  void logInit(const std::string& filename){
    logs::lout.open(filename.c_str());
    if(!logs::lout) {
      std::cerr << "error: log file cannot open." << std::endl;
      exit(1);
    }

    logs::lout << logs::info << "################ Program start ################" << std::endl;
  }

  void logConfig(const format::DockingConfiguration config){
    for(int i=0; i<config.ligand_files.size(); ++i){
      logs::lout << "Ligand file name    : "+config.ligand_files[i] << std::endl;
    }
    logs::lout << "Receptor file name  : "+config.receptor_file   << std::endl;
    logs::lout << "Output file name    : "+config.output_file     << std::endl;
    logs::lout << "Grid folder name    : "+config.grid_folder     << std::endl;
  }

  struct pos_param {
    int rotid, x, y, z;
    fltype score;
    int inp_ind;
    pos_param(int rotid, int x, int y, int z, fltype score, int inp_ind) : rotid(rotid), x(x), y(y), z(z), score(score), inp_ind(inp_ind) {}
    bool operator<(const pos_param& o) const { return score < o.score; }
    bool operator>(const pos_param& o) const { return score > o.score; }
  };

  // int to_search_num(int k, int score_num, int search_num, int ratio) {
  //   // assert(score_num >= search_num);
  //   // assert(0 <= k && k < score_num);
  //   // assert(0 <= search_num / 2 + (k - score_num / 2) / ratio && search_num / 2 + (k - score_num / 2) / ratio < search_num);
  //   return search_num / 2 + (k - score_num / 2) / ratio;
  // }

  int to_score_num(int k, int score_num, int search_num, int ratio) {
    // assert(score_num >= search_num);
    // assert(0 <= k && k < search_num);
    // assert(0 <= score_num / 2 + (k - search_num / 2) * ratio && score_num / 2 + (k - search_num / 2) * ratio < score_num);
    return score_num / 2 + (k - search_num / 2) * ratio;
  }

  fltype calcscore(const fragdock::Molecule& mol, const std::vector<fragdock::AtomEnergyGrid>& atom_grids) {
    fltype ret = 0.0;
    for (auto& a : mol.getAtoms()) {
      ret += atom_grids[a.getXSType()].getEnergy(a);
    }
    return ret;
  }

} // namespace

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;

  format::DockingConfiguration config = parseArgs(argc, argv);

  logInit(config.output_file + "fraggrid__" + getDate() + ".log");
  logConfig(config);

  // parse receptor file
  OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse ligands file
  vector<OpenBabel::OBMol> ligands = format::ParseFileToOBMol(config.ligand_files);

  int ligs_sz = ligands.size();
  logs::lout << "number of ligands: " << ligs_sz << endl;

  // ================================================================
  // prepare atomgrids and rotations
  // ================================================================
  logs::lout << logs::info << "[start] read energy grids" << endl;
  vector<AtomEnergyGrid> atom_grids = AtomEnergyGrid::readAtomGrids(config.grid_folder);
  logs::lout << logs::info << "[ end ] read energy grids" << endl;
  // logs::lout << "atom grid size: " << atom_grids.size() << endl;

  logs::lout << logs::info << "[start] make rotations" << endl;
  const vector<Vector3d> rotations = makeRotations60();
  logs::lout << logs::info << "[ end ] make rotations" << endl;


  const Point3d<int>& score_num = atom_grids[0].getNum();

  Point3d<int> ratio(utils::round(config.grid.search_pitch_x / config.grid.score_pitch_x),
                     utils::round(config.grid.search_pitch_y / config.grid.score_pitch_y),
                     utils::round(config.grid.search_pitch_z / config.grid.score_pitch_z));

  assert(abs(config.grid.score_pitch_x * ratio.x - config.grid.search_pitch_x) < 1e-4);
  assert(abs(config.grid.score_pitch_y * ratio.y - config.grid.search_pitch_y) < 1e-4);
  assert(abs(config.grid.score_pitch_z * ratio.z - config.grid.search_pitch_z) < 1e-4);

  Point3d<fltype> search_pitch(config.grid.search_pitch_x, config.grid.search_pitch_y, config.grid.search_pitch_z);
  Point3d<int> search_num(static_cast<int>(ceil(config.grid.inner_width_x / 2 / search_pitch.x) + 1e-9) * 2 + 1,
                          static_cast<int>(ceil(config.grid.inner_width_y / 2 / search_pitch.y) + 1e-9) * 2 + 1,
                          static_cast<int>(ceil(config.grid.inner_width_z / 2 / search_pitch.z) + 1e-9) * 2 + 1);


  EnergyGrid search_grid(atom_grids[0].getCenter(), search_pitch, search_num);

  vector<Molecule> ligands_mol(ligs_sz);

  int lig_kind_sz = 0;
  vector<vector<int> > same_ligs;
  unordered_map<string, int> lig_map;

  fltype margin = min({ config.grid.outer_width_x - config.grid.inner_width_x,
                        config.grid.outer_width_y - config.grid.inner_width_y,
                        config.grid.outer_width_z - config.grid.inner_width_z }) * 0.5;


  logs::lout << logs::info << "start pre-calculate energy" << endl;
  // for calc intra energy of ligands
  EnergyCalculator calc(1.0, 0.0);
  logs::lout << logs::info << "end pre-calculate energy" << endl;


  logs::lout << "grid margin : " << margin << endl;

  for (int lig_ind = 0; lig_ind < ligs_sz; ++lig_ind) {
    OpenBabel::OBMol& ligand = ligands[lig_ind];
    ligand.AddHydrogens();
    ligands_mol[lig_ind] = format::toFragmentMol(ligand);
    ligand.DeleteHydrogens();

    Molecule& mol = ligands_mol[lig_ind];

    mol.deleteHydrogens();
    mol.setIntraEnergy(calc.getIntraEnergy(mol));

    // mol.deleteHydrogens();
    // mol.calcRadius();

    Vector3d ofst = mol.getCenter();
    // logs::lout << ofst << endl;
    mol.translate(-ofst);

    const string& identifier = mol.getIdentifier();
    if (!lig_map.count(identifier)) {
      lig_map[identifier] = lig_kind_sz;
      ++lig_kind_sz;
      same_ligs.push_back(vector<int>());
    }
    same_ligs[lig_map[identifier]].push_back(lig_ind);
  }

  vector<int> sorted_lig(ligs_sz);
  for (int i = 0; i < ligs_sz; ++i) {
    sorted_lig[i] = i;
  }

  // vector<vector<OpenBabel::OBMol> > outputobmols(lig_kind_sz);

  int top_num = 1;
  int top_margin = 20;

  logs::lout << logs::info << top_num << " " << top_margin << endl;


  std::chrono::milliseconds fgrid_time(0);
  std::chrono::milliseconds real_time(0);

  vector<priority_queue<pos_param> > q(lig_kind_sz);

  int gsx = to_score_num(0, score_num.x, search_num.x, ratio.x);
  int gsy = to_score_num(0, score_num.y, search_num.y, ratio.y);
  int gsz = to_score_num(0, score_num.z, search_num.z, ratio.z);

  logs::lout << logs::info << "[TIME STAMP] START CALCULATING BY FRAGGRID" << endl;

  for (int i = 0; i < ligs_sz; ++i) {
    int lig_ind = sorted_lig[i];
    const Molecule& mol = ligands_mol[lig_ind];
    const string& identifier = mol.getIdentifier();
    // logs::lout << logs::info << (lig_ind + 1) << "th ligand : " << title << endl;

    int rotsz = rotations.size();

    vector<EnergyGrid> scores(rotsz, EnergyGrid(atom_grids[0].getCenter(), search_pitch, search_num, mol.getIntraEnergy()));


    auto t1 = std::chrono::system_clock::now();
    int loopcnt = 0;

      // #pragma omp parallel for
      for (int rotid = 0; rotid < rotsz; ++rotid) {
        Molecule mol = ligands_mol[lig_ind];
        mol.rotate(rotations[rotid]);

        for (int x = 0, gx = gsx; x < search_num.x; ++x, gx += ratio.x)
        for (int y = 0, gy = gsy; y < search_num.y; ++y, gy += ratio.y)
        for (int z = 0, gz = gsz; z < search_num.z; ++z, gz += ratio.z) {
          ++loopcnt;
          Molecule mmol = mol;
          mmol.translate(atom_grids[0].convert(gx, gy, gz));
          scores[rotid].addEnergy(x, y, z, calcscore(mmol, atom_grids));
        }
      }
    // logs::lout << logs::info << "end calc grid score" << endl;
    auto t2 = std::chrono::system_clock::now();

    fgrid_time += std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1);

    logs::lout << loopcnt << endl;

    // priority_queue<pos_param, vector<pos_param>, greater<pos_param> > q;
    int ind = lig_map[identifier];
    for (int rotid = 0; rotid < rotsz; ++rotid) {
      for (int x = 0; x < search_num.x; ++x) {
        for (int y = 0; y < search_num.y; ++y) {
          for (int z = 0; z < search_num.z; ++z) {
            q[ind].push(pos_param(rotid, x, y, z, scores[rotid].getEnergy(x, y, z), lig_ind));
            if (q[ind].size() > top_margin)
              q[ind].pop();
          }
        }
      }
    }

    auto t3 = std::chrono::system_clock::now();
    real_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t2);

    // logs::lout << logs::info << "end calc real score" << endl;
  }

  logs::lout << logs::info << "[TIME STAMP] END CALCULATING BY FRAGGRID" << endl;
  logs::lout << logs::info << "[TIME STAMP] fgrid_time : " << fgrid_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] realcalc_time : " << real_time.count() << endl;

  vector<pair<fltype, string>> ranking;
  ranking.reserve(lig_kind_sz);
  OpenBabel::outputOBMol outputs(config.output_file);
  ofstream outputcsv(config.output_file + "fraggrid__" + getDate() + ".csv");


  logs::lout << logs::info << "start pre-calculate energy" << endl;
  // for optimize
  calc = EnergyCalculator(1.0, -0.5);
  logs::lout << logs::info << "end pre-calculate energy" << endl;


  logs::lout << logs::info << "[TIME STAMP] START OPTIMIZING AND RANKING" << endl;

  for (const auto& p : lig_map) {
    const string& identifier = p.first;
    int i = p.second;
    vector<pair<fltype, pair<int, Molecule> > > out_mols;
    out_mols.reserve(q[i].size());
    while (!q[i].empty()) {
      pos_param param = q[i].top();
      q[i].pop();
      Molecule mol = ligands_mol[param.inp_ind];
      mol.rotate(rotations[param.rotid]);
      Vector3d pos = search_grid.convert(param.x, param.y, param.z);
      mol.translate(pos);
      Optimizer opt(receptor_mol);
      fltype opt_score = opt.optimize(mol, calc);

      out_mols.push_back(make_pair(opt_score, make_pair(param.inp_ind, mol)));
    }
    sort(out_mols.begin(), out_mols.end());

    fltype best_intra = ligands_mol[out_mols[0].second.first].getIntraEnergy();

    fltype best_score = out_mols[0].first / (1 + 0.05846 * ligands_mol[out_mols[0].second.first].getNrots());
    ranking.push_back(make_pair(best_score, identifier));

    outputcsv << identifier << "," << best_score << endl;

    // logs::lout << "mol : " << title << endl;

    for (int j = 0; j < top_num && j < out_mols.size(); ++j) {
      int lig_ind = out_mols[j].second.first;
      fltype score = (out_mols[j].first + ligands_mol[lig_ind].getIntraEnergy() - best_intra) / (1 + 0.05846 * ligands_mol[lig_ind].getNrots());
      logs::lout << "  " << (j + 1) << "th pose's score : " << score << endl;
      OpenBabel::OBMol obmol = ligands[lig_ind];
      OpenBabel::UpdateCoords(obmol, out_mols[j].second.second);
      outputs.write(obmol);

      // obmol.AddHydrogens();
      // Molecule mol = format::toFragmentMol(obmol);
    }
  }
  outputs.close();
  outputcsv.close();

  assert(ranking.size() == lig_kind_sz);
  sort(ranking.begin(), ranking.end());
  for (int i = 0; i < lig_kind_sz; ++i) {
    logs::lout << (i + 1) << "th ligand : " << ranking[i].second << endl;
    logs::lout << "score : " << ranking[i].first << endl;
  }

  logs::lout << logs::info << "[TIME STAMP] END OPTIMIZING AND RANKING" << endl;

  // OpenBabel::outputOBMolsToSDF(config.output_file, outputobmols);

  logs::lout << logs::info << "################ Program end ################" << endl;

  logs::lout << logs::info << "[FINAL_RESULT] "
      << "fgrid_time : " << fgrid_time.count() << ", "
      << "realcalc_time : " << real_time.count() << endl;

  // logs::lout << logs::info << "[RANKING_SCORE_RESULT]" << endl;
  // for (int i = 0; i < lig_kind_sz; ++i) {
  //   if (i > 0)
  //     logs::lout << " " << ranking[i].first;
  //   else
  //     logs::lout << ranking[i].first;
  // }
  // logs::lout << endl;

  logs::lout.close();
  return 0;
}
