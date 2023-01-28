#include "common.hpp"
#include "utils.hpp"
#include "OBMol.hpp"
#include "MoleculeToFragments.hpp"
#include "infile_reader.hpp"
#include "log_writer_stream.hpp"
#include "AtomEnergyGrid.hpp"
#include "FragmentEnergyGrid.hpp"
#include "FragmentsVector.hpp"
#include "EnergyCalculator.hpp"
#include "CalcMCFP.hpp"
#include "Optimizer.hpp"
#include "MinValuesVector.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <cstdlib>
#include <sys/stat.h>
#include <chrono>
#include <stdexcept>

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
      ("memsize,m", value<int64_t>(), "fragment grid's memory size[MB]")
      ("log", value<std::string>(), "log file");
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
    if (vmap.count("memsize")) conf.mem_size = vmap["memsize"].as<int64_t>();
    if (vmap.count("log")) conf.log_file = vmap["log"].as<std::string>();
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


  void logConfig(const format::DockingConfiguration config){
    for(int i=0; i<config.ligand_files.size(); ++i){
      logs::lout << "Ligand file name    : "+config.ligand_files[i] << std::endl;
    }
    logs::lout << "Receptor file name  : "+config.receptor_file   << std::endl;
    logs::lout << "Output file name    : "+config.output_file     << std::endl;
    logs::lout << "Grid folder name    : "+config.grid_folder     << std::endl;
    logs::lout << "Memory size[MB]     : "<<config.mem_size       << std::endl;
  }

  struct pos_param {
    int rotid, x, y, z;
    fltype score;
    int inp_ind;
    pos_param() : rotid(-1), x(-1), y(-1), z(-1), score(INF_ENERGY), inp_ind(-1) {}
    pos_param(int rotid, int x, int y, int z, fltype score, int inp_ind) : rotid(rotid), x(x), y(y), z(z), score(score), inp_ind(inp_ind) {}
    bool operator<(const pos_param& o) const { return score < o.score; }
    bool operator>(const pos_param& o) const { return score > o.score; }
  };

  std::vector<MCFP::node> makeGraph(const std::vector<fragdock::FragmentsVector>& frags, const std::vector<int>& sorted_lig, int fsize, int& reduce) {
    int sz = 0;
    for (auto& i : frags) {
      sz += i.size();
    }
    std::vector<int> before(fsize, -1);
    std::vector<MCFP::node> ret(sz, MCFP::node(-1, 0));
    int c = 0;
    for (int i = 0; i < frags.size(); ++i) {
      for (auto& j : frags[sorted_lig[i]].getvecs()) {
        int id = j.frag_idx;
        if (before[id] == c) {
          reduce += j.size;
          ret[c] = MCFP::node(c, 0);
        }
        else if (before[id] != -1)
          ret[before[id]] = MCFP::node(c, j.size);

        ++c;
        before[id] = c;
      }
    }
    return ret;
  }

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

  void inverse(std::vector<int>& v) {
    int n = v.size();
    std::vector<int> ret(n, -1);
    for (int i = 0; i < n; ++i) {
      assert(v[i] >= 0 && v[i] < n && ret[v[i]] == -1);
      ret[v[i]] = i;
    }
    for (int i = 0; i < n; ++i) {
      v[i] = ret[i];
    }
  }

} // namespace

namespace fragdock {
  EnergyGrid makeDistanceGrid(const Point3d<fltype>& center,
                              const Point3d<fltype>& pitch,
                              const Point3d<int>& num,
                              const Molecule& receptor_mol) {

    EnergyGrid grid(center, pitch, num);
    for(int x=0; x<grid.getNum().x; x++) {
      for(int y=0; y<grid.getNum().y; y++) {
        for(int z=0; z<grid.getNum().z; z++) {

          fltype mindist = 1e4;
          Vector3d pos = grid.convert(x, y, z);

          for (const Atom& a : receptor_mol.getAtoms()) {
            if (a.getXSType() != XS_TYPE_H) {
              mindist = std::min(mindist, (pos - a).abs());
            }
          }
          grid.setEnergy(x, y, z, mindist);
        }
      }
    }
    return grid;
  }

}

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;

  format::DockingConfiguration config = parseArgs(argc, argv);

  if(config.log_file == ""){
    config.log_file = config.output_file + "fraggrid__" + getDate() + ".log";
  }
  logs::log_init(config.log_file, true);
  logConfig(config);

  // parse receptor file
  OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse ligands file
  vector<OpenBabel::OBMol> ligands = format::ParseFileToOBMol(config.ligand_files);

  int ligs_sz = ligands.size(); /* the number of ligands */
  logs::lout << "number of ligands: " << ligs_sz << endl;

  // ================================================================
  // prepare atomgrids and rotations
  // ================================================================
  logs::lout << logs::info << "[start] read energy grids" << endl;
  vector<AtomEnergyGrid> atom_grids = AtomEnergyGrid::readAtomGrids(config.grid_folder);
  logs::lout << logs::info << "[ end ] read energy grids" << endl;
  // logs::lout << "atom grid size: " << atom_grids.size() << endl;

  vector<Vector3d> rotations_ligand;
  vector<Vector3d> rotations_fragment;
  if (config.rotangs_file.length()) {
    logs::lout << logs::info << "[start] read rotations" << endl;
    rotations_ligand = readRotations(config.rotangs_file);
    logs::lout << logs::info << "[ end ] read rotations" << endl;
    logs::lout << logs::info << "[start] make rotations" << endl;
    rotations_fragment = makeRotations60();
    logs::lout << logs::info << "[ end ] make rotations" << endl;
  }
  else {
    logs::lout << logs::info << "[start] make rotations" << endl;
    rotations_ligand = makeRotations60();
    rotations_fragment = rotations_ligand;
    logs::lout << logs::info << "[ end ] make rotations" << endl;
  }


  const Point3d<int>& score_num = atom_grids[0].getNum();

  Point3d<int> ratio(utils::round(config.grid.search_pitch_x / config.grid.score_pitch_x),
                     utils::round(config.grid.search_pitch_y / config.grid.score_pitch_y),
                     utils::round(config.grid.search_pitch_z / config.grid.score_pitch_z));

  // check whether the ratios are integer or not.
  assert(abs(config.grid.score_pitch_x * ratio.x - config.grid.search_pitch_x) < EPS);
  assert(abs(config.grid.score_pitch_y * ratio.y - config.grid.search_pitch_y) < EPS);
  assert(abs(config.grid.score_pitch_z * ratio.z - config.grid.search_pitch_z) < EPS);

  Point3d<fltype> search_pitch(config.grid.search_pitch_x, config.grid.search_pitch_y, config.grid.search_pitch_z);
  Point3d<int> search_num(utils::ceili(config.grid.inner_width_x / 2 / search_pitch.x) * 2 + 1,
                          utils::ceili(config.grid.inner_width_y / 2 / search_pitch.y) * 2 + 1,
                          utils::ceili(config.grid.inner_width_z / 2 / search_pitch.z) * 2 + 1);


  EnergyGrid search_grid(atom_grids[0].getCenter(), search_pitch, search_num);


  // The number how many fragment grids can be stored in memory.
  // It affects reuse of fragment grids.
  int FGRID_SIZE = (int)((config.mem_size * 1024 * 1024) / ((int64_t)score_num.x * score_num.y * score_num.z * sizeof(fltype)));
  if (FGRID_SIZE < 1) {
    logs::lout << logs::error << "Memory size is too small to load a fragment grid." << endl;
    logs::lout << logs::error << "Memory size of " << ((int64_t)score_num.x * score_num.y * score_num.z * sizeof(fltype) / 1024 / 1024) + 1 << " MB is at least needed." << endl;
    throw std::runtime_error("Memory size is too small to load a fragment grid.");
  }
  if (config.reuse_grid == format::DockingConfiguration::REUSE_NONE) {
    FGRID_SIZE = 1;
  }
  logs::lout << logs::info << "fragment grids storage size : " << FGRID_SIZE << endl;

  vector<FragmentsVector> fragvecs(ligs_sz); /* a vector of fragment vectors which correspond to ligands */

  vector<Fragment> frag_library; /* list of unique fragments which poses are normalized*/

  vector<int> frag_importance; /* # fragment heavy atoms * fragment occurrence */
  unordered_map<string, uint> fragmap; /* a map of {smiles -> an index of a fragment list} */

  vector<vector<uint> > same_ligs; /* conformer indices list for each unique ligand */
  unordered_map<string, uint> lig_map; /* smiles -> an index of a ligand list */

  logs::lout << logs::info << "start pre-calculate energy" << endl;
  EnergyCalculator ene_calculator(1.0); 
  logs::lout << logs::info << "end pre-calculate energy" << endl;


  logs::lout << logs::info << "[TIME STAMP] START MOLECULE OBJECT CONVERSION" << endl;
  vector<Molecule> ligands_mol(ligs_sz); /* a vector of ligand objects */
  for (uint l_ind = 0; l_ind < ligs_sz; ++l_ind) {
    OpenBabel::OBMol& ob_ligand = ligands[l_ind];
    ob_ligand.AddHydrogens();
    ligands_mol[l_ind] = format::toFragmentMol(ob_ligand);
    ob_ligand.DeleteHydrogens();
  }
  logs::lout << logs::info << "[TIME STAMP] END   MOLECULE OBJECT CONVERSION" << endl;


  logs::lout << logs::info << "[TIME STAMP] START INTRA ENERGY CALCULATION" << endl;
  for (uint l_ind = 0; l_ind < ligs_sz; ++l_ind) {
    Molecule& mol_ligand = ligands_mol[l_ind];
    mol_ligand.deleteHydrogens();
    mol_ligand.setIntraEnergy(
      ene_calculator.calcIntraEnergy(mol_ligand)
    );
  }
  logs::lout << logs::info << "[TIME STAMP] END   INTRA ENERGY CALCULATION" << endl;


  logs::lout << logs::info << "[TIME STAMP] START FRAGMENTATION" << endl;
  vector<vector<Fragment> > fragments_of_ligands(ligs_sz); /* list of fragments for each ligand */
  for (uint l_ind = 0; l_ind < ligs_sz; ++l_ind) {
    Molecule& mol_ligand = ligands_mol[l_ind];
    fragments_of_ligands[l_ind] = DecomposeMolecule(mol_ligand);
  }
  logs::lout << logs::info << "[TIME STAMP] END   FRAGMENTATION" << endl;




  for (uint l_ind = 0; l_ind < ligs_sz; ++l_ind) {
    OpenBabel::OBMol&  ob_ligand  = ligands[l_ind];
    Molecule&          mol_ligand = ligands_mol[l_ind];
    vector<Fragment>&  fragments  = fragments_of_ligands[l_ind];

    for (Fragment& mol_frag : fragments) {
      { // renumbering and setting smiles to mol_frag
        OpenBabel::OBMol ob_frag = format::toOBMol(mol_frag, ob_ligand);
        mol_frag.setSmiles(OpenBabel::canonicalSmiles(ob_frag));
        const vector<uint> canon_labels = OpenBabel::getRenumber(ob_frag); /* atom indices ordered canonically */
        mol_frag.renumbering(mol_frag.size(), canon_labels);
      } // destuct object ob_frag

      if (fragmap.count(mol_frag.getSmiles())) {
        uint frag_idx = fragmap[mol_frag.getSmiles()];
        mol_frag.setIdx(frag_idx);
        frag_importance[frag_idx] += mol_frag.size();
      } else {
        uint frag_idx = frag_library.size();
        mol_frag.setIdx(frag_idx);
        fragmap[mol_frag.getSmiles()] = frag_idx;

        Fragment temp = mol_frag;
        temp.normalize_pose();
        frag_library.push_back(temp);

        frag_importance.push_back(0);
      }

      fragvecs[l_ind].append(fragvec(mol_frag.getCenter(), mol_frag.getIdx(), mol_frag.size()));

      // logs::lout << temp << mol_frag;
    }

    Vector3d ofst = mol_ligand.getCenter();
    mol_ligand.translate(-ofst);
    fragvecs[l_ind].translate(-ofst);

    const string& identifier = mol_ligand.getIdentifier();
    if (!lig_map.count(identifier)) {
      lig_map[identifier] = lig_map.size();
      same_ligs.push_back(vector<uint>());
    }
    same_ligs[lig_map[identifier]].push_back(l_ind);
  }

  fragmap.clear();

  logs::lout << logs::info << "[TIME STAMP] END FRAGMENTATION" << endl;

  logs::lout << "config.reuse_grid : " << config.getReuseGridString() << endl;
  logs::lout << "config.reorder    : " << (config.reorder ? "True" : "False") << endl;


  vector<int> sorted_lig(ligs_sz);
  for (int i = 0; i < ligs_sz; ++i) {
    sorted_lig[i] = i;
  }

  // ---- use in REUSE_GRID_OFFLINE only ----
  vector<int> nextgridsp;
  int pred_reduce = 0;
  // ----------------------------------------

  logs::lout << logs::info << "[TIME STAMP] START REORDERING AND SOLVING MCFP" << endl;

  if (config.reorder) {
  // if (config.reuse_grid == format::DockingConfiguration::REUSE_OFFLINE) {
    vector<int> fragrank(frag_library.size());
    for (int i = 0; i < frag_library.size(); ++i) {
      fragrank[i] = i;
    }
    sort(fragrank.begin(), fragrank.end(), [&](const int& a, const int& b){ return frag_importance[a] > frag_importance[b]; });
    inverse(fragrank);
    // for (int i = 0; i < frag_library.size(); ++i) {
    //   logs::lout << fragrank[i] << " " << frag_importance[i] << endl;
    // }

    for (int i = 0; i < ligs_sz; ++i) {
      fragvecs[i].sort(fragrank);
    }
    sort(sorted_lig.begin(), sorted_lig.end(), [&](const int& a, const int& b){ return fragvecs[a] < fragvecs[b]; });
  }

  if (config.reuse_grid == format::DockingConfiguration::REUSE_OFFLINE) {
    vector<MCFP::node> graph = makeGraph(fragvecs, sorted_lig, frag_library.size(), pred_reduce);

    pred_reduce += MCFP::runLeftBackSSP(graph, FGRID_SIZE, nextgridsp);
    // cerr << "predict reduce cost : " << pred_reduce << endl;
  }

  logs::lout << logs::info << "[TIME STAMP] END REORDERING AND SOLVING MCFP" << endl;

  vector<FragmentEnergyGrid> fragment_grids(FGRID_SIZE);
  EnergyGrid distance_grid = makeDistanceGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol);

  // vector<vector<OpenBabel::OBMol> > outputobmols(lig_map.size());

  int top_num = 1;
  int top_before_strictopt = 2000;
  // int top_before_gridopt = 200;
  fltype SCORE_THRE = -3.0;

  int reduces = 0;

  // logs::lout << logs::info << top_num << " " << top_before_strictopt << ' ' << top_before_gridopt << endl;

  // ---- use in REUSE_GRID_ONLINE only ----
  vector<int> last_used(FGRID_SIZE, -1);
  vector<int> fgrid_ind(frag_library.size(), -1);
  // ---------------------------------------

  std::chrono::milliseconds fgrid_time(0);
  std::chrono::milliseconds real_time(0);

  // vector<priority_queue<pos_param> > q(lig_map.size());
  vector<utils::MinValuesVector<pos_param> > pos_param_vec(lig_map.size(), utils::MinValuesVector<pos_param>(top_before_strictopt));

  int gsx = to_score_num(0, score_num.x, search_num.x, ratio.x);
  int gsy = to_score_num(0, score_num.y, search_num.y, ratio.y);
  int gsz = to_score_num(0, score_num.z, search_num.z, ratio.z);

  logs::lout << logs::info << "[TIME STAMP] START CALCULATING BY FRAGGRID" << endl;

  for (int i = 0, frag_itr = 0; i < ligs_sz; ++i) {
    int lig_ind = sorted_lig[i];
    const Molecule& mol = ligands_mol[lig_ind];
    const string& identifier = mol.getIdentifier();
    // logs::lout << logs::info << (lig_ind + 1) << "th ligand : " << title << endl;

    int rotsz = rotations_ligand.size();
    int sz = fragvecs[lig_ind].size();

    vector<vector<int> > dx(rotsz, vector<int>(sz));
    vector<vector<int> > dy(rotsz, vector<int>(sz));
    vector<vector<int> > dz(rotsz, vector<int>(sz));

    vector<EnergyGrid> scores(rotsz, EnergyGrid(atom_grids[0].getCenter(), search_pitch, search_num, mol.getIntraEnergy()));
    // vector<EnergyGrid> scores(rotsz, EnergyGrid(atom_grids[0].getCenter(), search_pitch, search_num, 0.0));

    for (int rotid = 0; rotid < rotsz; ++rotid) {
      FragmentsVector fv = fragvecs[lig_ind];
      fv.rotate(rotations_ligand, rotid);

      for (int j = 0; j < sz; ++j) {
        dx[rotid][j] = utils::round(fv.getvec(j).pos.x / config.grid.score_pitch_x);
        dy[rotid][j] = utils::round(fv.getvec(j).pos.y / config.grid.score_pitch_y);
        dz[rotid][j] = utils::round(fv.getvec(j).pos.z / config.grid.score_pitch_z);
      }
    }

    auto t1 = std::chrono::system_clock::now();

    for (int j = 0; j < sz; ++j, ++frag_itr) {
      int fragid = fragvecs[lig_ind].getvec(j).frag_idx;
      int nextsp = -1;
      if (config.reuse_grid == format::DockingConfiguration::REUSE_OFFLINE || config.reuse_grid == format::DockingConfiguration::REUSE_NONE) {
        if (config.reuse_grid == format::DockingConfiguration::REUSE_OFFLINE)
          nextsp = nextgridsp[frag_itr];
        else
          nextsp = 0;

        if (fragment_grids[nextsp].frag_idx != fragid) {
          fragment_grids[nextsp] = FragmentEnergyGrid(frag_library[fragid], rotations_fragment, atom_grids, distance_grid);
        }
        else {
          reduces += fragvecs[lig_ind].getvec(j).size;
        }
      }
      else if (config.reuse_grid == format::DockingConfiguration::REUSE_ONLINE) {
        if (fgrid_ind[fragid] == -1) {
          int mi = frag_itr;
          for (int k = 0; k < FGRID_SIZE; ++k) {
            if (last_used[k] < mi) {
              mi = last_used[k];
              nextsp = k;
            }
          }
          assert(nextsp != -1);
          if (mi != -1) {
            // logs::lout << logs::info << "remove id : " << fragment_grids[nextsp].frag_idx << " fragment grid in vector at " << nextsp << endl;
            fgrid_ind[fragment_grids[nextsp].frag_idx] = -1;
          }
          // logs::lout << logs::info << "calc id : " << fragid << " fragment grid and store to vector at " << nextsp << endl;
          fragment_grids[nextsp] = FragmentEnergyGrid(frag_library[fragid], rotations_fragment, atom_grids, distance_grid);
          fgrid_ind[fragid] = nextsp;
        }
        else {
          nextsp = fgrid_ind[fragid];
          // logs::lout << logs::info << "get id : " << fragid << " fragment grid from vector at " << nextsp << endl;
          reduces += fragvecs[lig_ind].getvec(j).size;
        }
        last_used[fgrid_ind[fragid]] = frag_itr;
      }
      else {
        assert(0);
      }

      assert(fragment_grids[nextsp].frag_idx == fragid);
      const FragmentEnergyGrid& fg = fragment_grids[nextsp];

      // #pragma omp parallel for
      for (int rotid = 0; rotid < rotsz; ++rotid) {
        // int rid = RotMatrix[fragvecs[lig_ind].getvec(j).rotid][rotid];

        for (int x = 0, gx = gsx; x < search_num.x; ++x, gx += ratio.x)
        for (int y = 0, gy = gsy; y < search_num.y; ++y, gy += ratio.y)
        for (int z = 0, gz = gsz; z < search_num.z; ++z, gz += ratio.z)
          scores[rotid].addEnergy(x, y, z, fg.getGrid().getEnergy(gx + dx[rotid][j], gy + dy[rotid][j], gz + dz[rotid][j]));
      }
    }
    // logs::lout << logs::info << "end calc grid score" << endl;
    auto t2 = std::chrono::system_clock::now();

    fgrid_time += std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1);

    // priority_queue<pos_param, vector<pos_param>, greater<pos_param> > q;
    int ind = lig_map[identifier];
    for (int rotid = 0; rotid < rotsz; ++rotid) {
      for (int x = 0; x < search_num.x; ++x) {
        for (int y = 0; y < search_num.y; ++y) {
          for (int z = 0; z < search_num.z; ++z) {
            const fltype score = scores[rotid].getEnergy(x, y, z);
            if (score < SCORE_THRE) {
              pos_param_vec[ind].push(pos_param(rotid, x, y, z, score, lig_ind));
              // q[ind].push(pos_param(rotid, x, y, z, score, lig_ind));
              // if (q[ind].size() > top_before_strictopt)
              //   q[ind].pop();
            }
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
  ranking.reserve(lig_map.size());
  OpenBabel::outputOBMol outputs(config.output_file);
  ofstream outputcsv(config.output_file + "fraggrid__" + getDate() + ".csv");


  logs::lout << logs::info << "start pre-calculate energy" << endl;
  // for optimize
  // calc = EnergyCalculator(0.95, 0.0);

  // Optimizer opt(receptor_mol);
  // vector<AtomEnergyGrid> opt_atom_grids = AtomEnergyGrid::makeAtomGrids(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol, calc);
  Optimizer_Grid opt_grid(atom_grids);

  logs::lout << logs::info << "end pre-calculate energy" << endl;


  logs::lout << logs::info << "[TIME STAMP] START OPTIMIZING AND RANKING" << endl;

  for (const auto& p : lig_map) {
    const string& identifier = p.first;
    int i = p.second;
    const vector<pos_param>& poss = pos_param_vec[i].getValues();
    vector<pair<fltype, pair<int, Molecule> > > out_mols(poss.size());
    for (int j = 0; j < (int)poss.size(); j++) {
      const pos_param& param = poss[j];
      Molecule mol = ligands_mol[param.inp_ind];
      mol.rotate(rotations_ligand[param.rotid]);
      Vector3d pos = search_grid.convert(param.x, param.y, param.z);
      mol.translate(pos);
      fltype opt_score = opt_grid.optimize(mol);
      // fltype opt_score = param.score;

      out_mols[j] = make_pair(opt_score, make_pair(param.inp_ind, mol));
    }
    sort(out_mols.begin(), out_mols.end());
    fltype best_intra = LIMIT_ENERGY;
    fltype best_score = LIMIT_ENERGY;
    if (out_mols.size() > 0) {
      best_intra = ligands_mol[out_mols[0].second.first].getIntraEnergy();
      best_score = out_mols[0].first / (1 + 0.05846 * ligands_mol[out_mols[0].second.first].getNrots());
    }
    ranking.push_back(make_pair(best_score, identifier));

    outputcsv << identifier << "," << best_score << endl;

    // logs::lout << "mol : " << title << endl;

    for (int j = 0; j < top_num && j < out_mols.size(); ++j) {
      int lig_ind = out_mols[j].second.first;
      fltype score = (out_mols[j].first + ligands_mol[lig_ind].getIntraEnergy() - best_intra) / (1 + 0.05846 * ligands_mol[lig_ind].getNrots());
      // fltype score = (out_mols[j].first) / (1 + 0.05846 * ligands_mol[lig_ind].getNrots());
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

  assert(ranking.size() == lig_map.size());
  sort(ranking.begin(), ranking.end());
  for (int i = 0; i < lig_map.size(); ++i) {
    logs::lout << (i + 1) << "th ligand : " << ranking[i].second << endl;
    logs::lout << "score : " << ranking[i].first << endl;
  }

  logs::lout << logs::info << "[TIME STAMP] END OPTIMIZING AND RANKING" << endl;

  // logs::lout << "predict reduce cost : " << pred_reduce << endl;
  logs::lout << "real reduce cost    : " << reduces << endl;
  // cerr << "real reduce cost    : " << reduces << endl;
  assert(config.reuse_grid != format::DockingConfiguration::REUSE_OFFLINE || reduces == pred_reduce);

  // OpenBabel::outputOBMolsToSDF(config.output_file, outputobmols);

  logs::lout << logs::info << "################ Program end ################" << endl;


  /* ########### calculate statistics ################ */
  long long allcost = 0; /* The total cost without the reuse of fragment grids */
  for (const auto& fragments: fragments_of_ligands) {
    for (const Fragment& mol_frag: fragments) {
      allcost += mol_frag.size();
    }
  }

  long long minimum_cost = 0; /* The total cost if all fragment grids can be stored */
  for (const auto& frag: frag_library){
    minimum_cost += frag.size();
  } 

  int fnums = 0;
  for (const auto& fragments: fragments_of_ligands) {
    fnums += fragments.size();
  }
  /* ########### calculate statistics END ################ */




  logs::lout << logs::info << "[FINAL_RESULT] "
      << ligs_sz << "/" << config.getReuseGridString() << "/" << (config.reorder ? "reorder" : "no_reorder") << "/" << config.mem_size << "/" << config.grid.inner_width_x << ", "
      << "fragment types : " << frag_library.size() << ", "
      << "fragment num : " << fnums << ", "
      << "all cost : " << allcost << ", "
      << "minimum cost : " << minimum_cost << ", "
      << "reduce cost : " << reduces << ", "
      << "real cost : " << allcost - reduces << ", "
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

  logs::close();
  return 0;
}
