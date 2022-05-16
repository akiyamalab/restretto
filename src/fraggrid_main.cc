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
#include "GradientCalculator.hpp"
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
#include <queue>

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
      ("memsize,m", value<int64_t>(), "fragment grid's memory size[MB]")//;
      
      ("maxiterations", value<int64_t>(), "max iterations")
      ("stepsize", value<fltype>(), "step size");
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
    
    if (vmap.count("maxiterations")) conf.max_iterations = vmap["maxiterations"].as<int64_t>();
    if (vmap.count("stepsize")) conf.step_size = vmap["stepsize"].as<fltype>();
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
      std::cout << std::endl << "Ligand file name          : "+config.ligand_files[i] << std::endl;
    }
    logs::lout << "Receptor file name        : "+config.receptor_file   << std::endl;
    logs::lout << "Output file name          : "+config.output_file     << std::endl;
    logs::lout << "Grid folder name          : "+config.grid_folder     << std::endl;
    logs::lout << "Fragment Grid folder name : "+config.fragment_grid_folder     << std::endl;
    logs::lout << "Memory size[MB]           : "<<config.mem_size       << std::endl;
    
    logs::lout << "Max iterations            : "<<config.max_iterations << std::endl;
    logs::lout << "Step size                 : "<<config.step_size      << std::endl;
    
    std::cout << "Receptor file name        : "+config.receptor_file   << std::endl;
    std::cout << "Output file name          : "+config.output_file     << std::endl;
    std::cout << "Grid folder name          : "+config.grid_folder     << std::endl;
    std::cout << "Fragment Grid folder name : "+config.fragment_grid_folder     << std::endl;
    std::cout << "Memory size[MB]           : "<<config.mem_size       << std::endl;
    
    std::cout << "Max iterations            : "<<config.max_iterations << std::endl;
    std::cout << "Step size                 : "<<config.step_size      << std::endl;
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
        int id = j.fragid;
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



  std::vector<std::pair<int, fltype>> makeFgridsBestScore(const std::vector<FragmentEnergyGrid>& fgrids) {

    int num_frags = fgrids.size();
    std::vector<std::pair<int, fltype>> fgrids_best_score(num_frags);
    for(int i = 0; i < num_frags; i++) {
      fgrids_best_score[i] = std::make_pair(i, fgrids[i].getBestScore());
    }
    std::sort(fgrids_best_score.begin(), fgrids_best_score.end(), [&](const std::pair<int, fltype>& p, const std::pair<int, fltype>& q){ return p.second < q.second; });

    return fgrids_best_score;
  }
  
  
  std::vector<int> frag_selection(const std::vector<std::pair<int, fltype>>& fgrids_best_score, const int LIM_SELE_FRAG_SZ, const fltype th_f) {

    std::vector<int> sele_frags;
    for(int i = 0; i < LIM_SELE_FRAG_SZ; i++) {
      if(fgrids_best_score[i].second >= th_f)
        break;
      sele_frags.push_back(fgrids_best_score[i].first);
    }

    return sele_frags;
  }
  
  
  std::vector<std::pair<int, fltype>> makeF1(const std::vector<std::pair<int, fltype>>& fgrids_best_score, const int SELE_FRAG_SZ) {

    std::vector<std::pair<int, fltype>> f1;
    for(int i = 0; i < SELE_FRAG_SZ; i++) {
      f1.push_back(fgrids_best_score[i]);
    }

    return f1;
  }


  std::vector<std::pair<int, fltype>> lig_selection(const std::vector<std::pair<int, fltype>>& f1, const std::vector<std::vector<std::pair<int, int>>>& frag_lig_list, const std::vector<int>& orig_ligs, const int LIM_SELE_LIG_SZ, const fltype th_l) {
    
    int num_ligs = orig_ligs.size();
    std::vector<std::pair<int, fltype>> ligs_score(num_ligs);
    std::map<int, int> dict;
    for (int i = 0; i < num_ligs; ++i) {
      ligs_score[i] = std::make_pair(orig_ligs[i], LIMIT_ENERGY);
      assert(dict.count(orig_ligs[i])==0);
      dict[orig_ligs[i]] = i;
    }

    for(int frag_itr = 0; frag_itr < f1.size(); frag_itr++) {
      int fragid = f1[frag_itr].first;
      fltype f1_score = f1[frag_itr].second; // not fragid but frag_itr

      for(int i = 0; i < frag_lig_list[fragid].size(); i++) {
        int lig_ind = frag_lig_list[fragid][i].first;

        assert(ligs_score[dict[lig_ind]].first == lig_ind);
        if(ligs_score[dict[lig_ind]].second > f1_score) {
          ligs_score[dict[lig_ind]].second = f1_score;
        }
      }
    }
    std::sort(ligs_score.begin(), ligs_score.end(), [&](const std::pair<int, fltype>& p, const std::pair<int, fltype>& q){ return p.second < q.second; });


    std::vector<std::pair<int, fltype>> sele_ligs_score;
    for(int i = 0; i < LIM_SELE_LIG_SZ; i++) {
      if(ligs_score[i].second >= th_l)
        break;
      sele_ligs_score.push_back(ligs_score[i]);
    }
    
    return sele_ligs_score;
  }



  std::map<std::pair<int, int>, std::vector<fltype>> makeF2(const std::vector<FragmentEnergyGrid>& fgrids, const std::set<std::pair<int, int>>& f2pair, const int top_fgrid_num, const int DMAX) {
    
    std::map<std::pair<int, int>, std::vector<fltype>> f2;

    for(auto itr = f2pair.begin(); itr != f2pair.end(); itr++) {
      f2[*itr] = std::vector<fltype>(DMAX, LIMIT_ENERGY);
      int fa_id = itr->first;
      int fb_id = itr->second;
      assert(fa_id <= fb_id);
      
      const std::vector<XYZ_param>& xyz_a = fgrids[fa_id].getGoodXYZs();
      assert(top_fgrid_num <= xyz_a.size());
      const std::vector<XYZ_param>& xyz_b = fgrids[fb_id].getGoodXYZs();
      assert(top_fgrid_num <= xyz_b.size());
      
      for(int s = 0; s < top_fgrid_num; s++) {
        const XYZ_param& pa = xyz_a[s];
        for(int t = 0; t < top_fgrid_num; t++) {
          const XYZ_param& pb = xyz_b[t];
            
          int dist = utils::round((pa.pos-pb.pos).abs());
          assert(dist < DMAX);

          if(f2[std::make_pair(fa_id, fb_id)][dist] > pa.score + pb.score) {
            f2[std::make_pair(fa_id, fb_id)][dist] = pa.score + pb.score;
          } 
        }
      }
        
    }

    return f2;
  }



  std::vector<std::pair<int, fltype>> conf_selection(const std::vector<std::pair<int, fltype>>& f1, const std::map<std::pair<int, int>, std::vector<fltype>>& f2, const std::vector<FragmentsVector>& fragvecs, const int DMAX, const fltype score_pitch, const std::vector<int>& orig_confs, const int LIM_SELE_CONF_SZ, const fltype th_c) {

    int num_confs = orig_confs.size();
    std::vector<std::pair<int, fltype>> confs_score(num_confs);
    for (int i = 0; i < num_confs; ++i) {
      confs_score[i] = std::make_pair(orig_confs[i], LIMIT_ENERGY);
    }

    std::map<int, int> revf1;
    for(int i = 0; i < f1.size(); i++) {
      revf1[f1[i].first] = i;
    }

    for (int i = 0; i < num_confs; ++i) {
      int conf_ind = confs_score[i].first;

      fltype best_f2_score = LIMIT_ENERGY;

      const FragmentsVector& fv = fragvecs[conf_ind];
      int sz = fv.size();
      
      if(sz == 1) {
        int existInF1 = revf1.count(fv.getvec(0).fragid);
        assert(existInF1);
        confs_score[i].second = f1[revf1[fv.getvec(0).fragid]].second;
      
      }else {
        for(int j = 0; j < sz; j++) {
          int fj_id = fv.getvec(j).fragid;
          Vector3d vj = fv.getvec(j).pos;

          for(int k = j + 1; k < sz; k++) {
            int fk_id = fv.getvec(k).fragid;
            Vector3d vk = fv.getvec(k).pos;

            int id0 = std::min(fj_id, fk_id);
            int id1 = std::max(fj_id, fk_id);
            int dist = utils::round((vj-vk).abs() / score_pitch);

            if(dist >= DMAX) continue;
            
            auto it = f2.find(std::make_pair(id0, id1));
            if(it == f2.end()) continue;
            
            fltype f2_score = it->second[dist];
            if(best_f2_score > f2_score) {
              best_f2_score = f2_score;
            }
          }
        }
      
        confs_score[i].second = best_f2_score;
      }

    }
    std::sort(confs_score.begin(), confs_score.end(), [&](const std::pair<int, fltype>& p, const std::pair<int, fltype>& q){ return p.second < q.second; });

    std::vector<std::pair<int, fltype>> sele_confs_score;
    for(int i = 0; i < LIM_SELE_CONF_SZ; i++) {
      if(confs_score[i].second >= th_c)
        break;
      sele_confs_score.push_back(confs_score[i]);
    }
    
    return sele_confs_score;
  }



}

namespace {

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

  using namespace std;
  void printDatasetInfo(const vector<fragdock::Molecule>& ligands_mol, const vector<fragdock::Fragment>& fragtemps, const vector<fragdock::FragmentsVector>& fragvecs, const vector<vector<int>>& same_ligs, const vector<vector<pair<int, int>>>& frag_lig_list, const vector<pair<int, fltype>>& fgrids_best_score, const vector<string>& fragsmiles, const int lig_kind_sz, const int ligs_sz, const int ftmpsz, const int fnums) {
    cout << "fgrids_best_score ranking:";
    for(int i = 0; i < fgrids_best_score.size(); i++) {
      cout << fgrids_best_score[i].first << ",";
    }
    cout << endl;
  
    cout << "fgrids_best_score score:";
    for(int i = 0; i < fgrids_best_score.size(); i++) {
      cout << fgrids_best_score[i].second << ",";
    }
    cout << endl;

    // 上位 rate % のフラグメントを1つ，2つ持っている化合物の割合
    cout << "express ligand:" << endl;
    vector<fltype> x;
    vector<fltype> y1, y2;
    for(int rate = 0; rate <= 100; rate++) {
      vector<int> cnt(lig_kind_sz, 0);
      for(int i = 0; i < utils::round(fgrids_best_score.size() * rate / 100.0); i++) {
        int fragid = fgrids_best_score[i].first;
    
        for(int j = 0; j < frag_lig_list[fragid].size(); j++) {
          cnt[frag_lig_list[fragid][j].first]++;
        }
      }
  
      int c1=0;
      int c2=0;
      for(int i = 0; i < lig_kind_sz; i++) {
        if(cnt[i]>=1) c1++;
        if(cnt[i]>=2) c2++;
      }
    
      x.push_back(rate);
      y1.push_back(fltype(c1)/lig_kind_sz * 100);
      y2.push_back(fltype(c2)/lig_kind_sz * 100);
    }
  
    for(auto v:  x) cout << v << ","; cout << endl;
    for(auto v: y1) cout << v << ","; cout << endl;
    for(auto v: y2) cout << v << ","; cout << endl;
  
  
    cout << endl;
    // 化合物のID, 持っているフラグメント数, タイトル，SMILES
    cout << "//ligid,atomnum,fnum,frags,ligtitle,ligsmiles" << endl;
    for(int i = 0; i < lig_kind_sz; i++) {
      cout << "ligid " << i << "," << ligands_mol[same_ligs[i][0]].size() << "," << fragvecs[same_ligs[i][0]].size() << ",";
    
      for(int j = 0; j < fragvecs[same_ligs[i][0]].size(); j++) {
        cout << fragvecs[same_ligs[i][0]].getvec(j).fragid << ",";
      }
      cout << ligands_mol[same_ligs[i][0]].getIdentifier() << endl;
    }
    cout << endl;

    // フラグメントID, 原子数, SMILES
    cout << "//fragid,atomnum,fragsmiles" << endl;
    for(int i = 0; i < ftmpsz; i++) {
      cout << "fragid " << i << "," << fragtemps[i].size() << "," << fragsmiles[i] << endl;
    }
    cout << endl << endl;
    
    //cout << "ligtypes,conftypes,fragtypes,fnums: " << lig_kind_sz << "," << ligs_sz << "," << ftmpsz << "," << fnums << endl;
  }

}

int main(int argc, char **argv){
  using namespace std;
  using namespace fragdock;

  format::DockingConfiguration config = parseArgs(argc, argv);

  //logInit(config.output_file + "fraggrid__" + getDate() + ".log");
  logInit(config.output_file + "_fraggrid" + ".log");
  logConfig(config);

  // parse receptor file
  OpenBabel::OBMol receptor = format::ParseFileToOBMol(config.receptor_file.c_str())[0];
  Molecule receptor_mol = format::toFragmentMol(receptor);

  // parse ligands file
  vector<OpenBabel::OBMol> ligands = format::ParseFileToOBMol(config.ligand_files);

  int ligs_sz = ligands.size();
  //logs::lout << "number of ligands: " << ligs_sz << endl;

  // ================================================================
  // prepare atomgrids and rotations
  // ================================================================
  logs::lout << logs::info << "[start] read energy grids" << endl;
  vector<AtomEnergyGrid> atom_grids = AtomEnergyGrid::readAtomGrids(config.grid_folder);
  logs::lout << logs::info << "[ end ] read energy grids" << endl;
  // logs::lout << "atom grid size: " << atom_grids.size() << endl;

  vector<vector<AtomEnergyGrid> > gradient_grids = GradientCalculator::readGradientGrids(config.grid_folder);
  assert(config.fragment_grid_folder.size()); // fragment_grid_folder

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

  assert(abs(config.grid.score_pitch_x * ratio.x - config.grid.search_pitch_x) < 1e-4);
  assert(abs(config.grid.score_pitch_y * ratio.y - config.grid.search_pitch_y) < 1e-4);
  assert(abs(config.grid.score_pitch_z * ratio.z - config.grid.search_pitch_z) < 1e-4);

  Point3d<fltype> search_pitch(config.grid.search_pitch_x, config.grid.search_pitch_y, config.grid.search_pitch_z);
  Point3d<int> search_num(static_cast<int>(ceil(config.grid.inner_width_x / 2 / search_pitch.x) + 1e-9) * 2 + 1,
                          static_cast<int>(ceil(config.grid.inner_width_y / 2 / search_pitch.y) + 1e-9) * 2 + 1,
                          static_cast<int>(ceil(config.grid.inner_width_z / 2 / search_pitch.z) + 1e-9) * 2 + 1);


  EnergyGrid search_grid(atom_grids[0].getCenter(), search_pitch, search_num);


  int FGRID_SIZE = (int)((config.mem_size << 20) / ((int64_t)score_num.x * score_num.y * score_num.z * sizeof(fltype)));
  if (FGRID_SIZE < 1) {
    cerr << "[assert] memory size is so small." << endl;
    cerr << "         need memory size : " << ((int64_t)score_num.x * score_num.y * score_num.z * sizeof(fltype) >> 20) + 1 << " MB" << endl;
    assert(0);
  }
  if (config.reuse_grid == format::DockingConfiguration::REUSE_NONE) {
    FGRID_SIZE = 1;
  }
  logs::lout << logs::info << "fragment grids storage size : " << FGRID_SIZE << endl;

  vector<Molecule> ligands_mol(ligs_sz);
  vector<FragmentsVector> fragvecs(ligs_sz);

  vector<Fragment> fragtemps;
  vector<string> fragsmiles; //新たに追加

  vector<int> frag_importance;
  int ftmpsz = 0;
  int fnums = 0;
  unordered_map<string, int> fragmap;

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


  logs::lout << logs::info << "[TIME STAMP] START FRAGMENTATION" << endl;

  logs::lout << "grid margin : " << margin << endl;
  long long allcost = 0;
  long long minimum_cost = 0;
  for (int lig_ind = 0; lig_ind < ligs_sz; ++lig_ind) { //フラグメント分割
    OpenBabel::OBMol& ligand = ligands[lig_ind];
    ligand.AddHydrogens();
    ligands_mol[lig_ind] = format::toFragmentMol(ligand);
    ligand.DeleteHydrogens();

    Molecule& mol = ligands_mol[lig_ind];

    mol.deleteHydrogens();
    fltype intra = calc.getIntraEnergy(mol);
    mol.setIntraEnergy(intra);

    vector<Fragment> fragments = DecomposeMolecule(mol);

    // logs::lout << (lig_ind + 1) << "th ligand, fragment size : " << fragments.size() << endl;

    for (auto& frag : fragments) {
      ++fnums;
      OpenBabel::OBMol obmol = format::toOBMol(frag, ligand);
      vector<unsigned int> canon_labels;
      OpenBabel::getRenumber(obmol, canon_labels);

      string smiles = OpenBabel::canonicalSmiles(obmol);

      frag.renumbering(frag.size(), canon_labels);
      allcost += frag.size();
      Fragment temp = frag;

      if (fragmap.count(smiles)) { // 同じフラグメントが他の化合物で既出
        temp = fragtemps[fragmap[smiles]];
        frag_importance[fragmap[smiles]] += temp.size();
      }
      else {
        temp.settri();
        temp.normalize();
        temp.settempind(ftmpsz);
        fragtemps.push_back(temp);
        
        fragsmiles.push_back(smiles);
        
        frag_importance.push_back(0);
        fragmap[smiles] = ftmpsz; //そのフラグメントのidの役割
        ++ftmpsz;
        minimum_cost += temp.size();
      }
      frag.settri(temp);
      // frag.setrotid(rotations);
      frag.settempind(temp.gettempind());

      fragvecs[lig_ind].append(fragvec(frag.getCenter(), frag.gettempind(), frag.size()));

      // logs::lout << temp << frag;
    }

    // mol.deleteHydrogens();
    // mol.calcRadius();
    Vector3d ofst = mol.getCenter();
    // logs::lout << ofst << endl;
    mol.translate(-ofst); //化合物の重心が原点になるように移動(*)
    fragvecs[lig_ind].translate(-ofst); //各フラグメントを(*)に合わせて移動

    const string& identifier = mol.getIdentifier(); //同じ化合物から出来た配座は同じidentifierを持つ
    if (!lig_map.count(identifier)) {
      lig_map[identifier] = lig_kind_sz;
      ++lig_kind_sz;
      same_ligs.push_back(vector<int>());
    }
    same_ligs[lig_map[identifier]].push_back(lig_ind); //identifierごとに配座のインデックスをまとめたもの．
  }

  //logs::lout << logs::info << "fragment types : " << ftmpsz << endl;

  //logs::lout << logs::info << "fragment num : " << fnums << endl;


  fragmap.clear();

  logs::lout << logs::info << "[TIME STAMP] END FRAGMENTATION" << endl;

  logs::lout << "fragment num         : " << fnums << endl;
  logs::lout << "fragment types       : " << ftmpsz << endl;
  logs::lout << "number of ligands    : " << lig_kind_sz << endl;
  logs::lout << "number of conformers : " << ligs_sz << endl;
  logs::lout << endl;

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
    vector<int> fragrank(ftmpsz);
    for (int i = 0; i < ftmpsz; ++i) {
      fragrank[i] = i;
    }
    sort(fragrank.begin(), fragrank.end(), [&](const int& a, const int& b){ return frag_importance[a] > frag_importance[b]; });
    inverse(fragrank);
    // for (int i = 0; i < ftmpsz; ++i) {
    //   logs::lout << fragrank[i] << " " << frag_importance[i] << endl;
    // }

    for (int i = 0; i < ligs_sz; ++i) {
      fragvecs[i].sort(fragrank); //化合物内のフラグメントを重要度順に並び替える
    }
    sort(sorted_lig.begin(), sorted_lig.end(), [&](const int& a, const int& b){ return fragvecs[a] < fragvecs[b]; }); //重要度の高いフラグメント集合を持っている化合物ほど前にくるようにソート
  }

  if (config.reuse_grid == format::DockingConfiguration::REUSE_OFFLINE) {
    vector<MCFP::node> graph = makeGraph(fragvecs, sorted_lig, ftmpsz, pred_reduce);

    pred_reduce += MCFP::runLeftBackSSP(graph, FGRID_SIZE, nextgridsp);
    cerr << "predict reduce cost : " << pred_reduce << endl;
  }

  logs::lout << logs::info << "[TIME STAMP] END REORDERING AND SOLVING MCFP" << endl;

  // vector<FragmentEnergyGrid> fragment_grids(FGRID_SIZE);
  vector<FragmentEnergyGrid> fragment_grids(ftmpsz);
  EnergyGrid distance_grid = makeDistanceGrid(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol);

  // vector<vector<OpenBabel::OBMol> > outputobmols(lig_kind_sz);

  int top_num = 1;
  int top_before_strictopt = 100;
  // int top_before_gridopt = 200;
  //fltype SCORE_THRE = -3.0;        //閾値
  fltype SCORE_THRE = 0.0;        //閾値
  int reduces = 0;

  // logs::lout << logs::info << top_num << " " << top_before_strictopt << ' ' << top_before_gridopt << endl;

  // ---- use in REUSE_GRID_ONLINE only ----
  //vector<int> last_used(FGRID_SIZE, -1);
  //vector<int> fgrid_ind(ftmpsz, -1);
  // ---------------------------------------

  //std::chrono::milliseconds fgrid_time(0);
  //std::chrono::milliseconds real_time(0);
  std::chrono::milliseconds fgrid_time(0);

  // vector<priority_queue<pos_param> > q(lig_kind_sz);
  vector<utils::MinValuesVector<pos_param> > pos_param_vec(lig_kind_sz, utils::MinValuesVector<pos_param>(top_before_strictopt));

  int gsx = to_score_num(0, score_num.x, search_num.x, ratio.x);
  int gsy = to_score_num(0, score_num.y, search_num.y, ratio.y);
  int gsz = to_score_num(0, score_num.z, search_num.z, ratio.z);



//-----------------------
  
  logs::lout << logs::info << "[TIME STAMP] START REORDERING BY LIGAND LEVEL" << endl;

  // 化合物レベルで reorder
  vector<int> original_ligs(lig_kind_sz); // ligs_sz ではない
  for (int i = 0; i < lig_kind_sz; ++i) { // ligs_sz ではない
    original_ligs[i] = i;
  }

  if (config.reorder) {
    sort(original_ligs.begin(), original_ligs.end(), [&](const int& a, const int& b){ return fragvecs[same_ligs[a][0]] < fragvecs[same_ligs[b][0]]; }); // 化合物レベルでのソートをするためsame_ligs[a], same_ligs[b] のうちそれぞれ代表の配座一つのfragvecを取ってきて比べる．
  }

  logs::lout << logs::info << "[TIME STAMP] END REORDERING BY LIGAND LEVEL" << endl;
//-----------------------
  logs::lout << logs::info << "[TIME STAMP] START MAKE FRAG_LIG_LIST" << endl;

  // 0. frag_lig_list 生成 
  vector<vector<pair<int, int>> > frag_lig_list(ftmpsz, vector<pair<int, int>>());
  for(int i = 0; i < lig_kind_sz; i++) {
    int conf_ind = same_ligs[original_ligs[i]][0];

    FragmentsVector& fv = fragvecs[conf_ind];
    int sz = fv.size();
    for(int j = 0; j < sz; j++) {
      frag_lig_list[fv.getvec(j).fragid].push_back(make_pair(original_ligs[i], j));
    }
  }

  logs::lout << logs::info << "[TIME STAMP] END MAKE FRAG_LIG_LIST" << endl;
//-----------------------
  logs::lout << logs::info << "[TIME STAMP] START MAKE FRAGGRID" << endl;

  // 1. フラグメントグリッド生成
  auto s1 = std::chrono::system_clock::now();
  int MAXK = 100;
  
  
  map<char, string> repla_dict;
  repla_dict['*'] = "a_";
  repla_dict['/'] = "s_";
  repla_dict['\\'] = "b_";
  repla_dict['#'] = "h_";
  
  repla_dict['+'] = "p_";
  repla_dict['-'] = "m_";
  repla_dict['='] = "e_";
  
  repla_dict['['] = "lsb_";
  repla_dict[']'] = "rsb_";
  
  repla_dict['('] = "lpa_";
  repla_dict[')'] = "rpa_";
  
  repla_dict['@'] = "at_";
  
  makeFolder(config.fragment_grid_folder);
  for(int i = 0; i < ftmpsz; i++) {
    
    cerr << fragsmiles[i] << endl;
    
    string replaced_str = "";
    for(auto c: fragsmiles[i]) {
      
      if(repla_dict.count(c)) {
        replaced_str += repla_dict[c];
      }else {
        replaced_str += c;
      }
      
    }
    
    string fragment_grid_filename = config.fragment_grid_folder + "/" + replaced_str + ".grid";
    ifstream ifs(fragment_grid_filename.c_str(), ios::binary);
    
    if(replaced_str.size() > 245) {
      cerr << "smiles too long. so not save fragmentgrid file " << fragsmiles[i] << endl;
      
      fragment_grids[i] = FragmentEnergyGrid(fragtemps[i], rotations_fragment, atom_grids, distance_grid, MAXK);
      assert(fragment_grids[i].getGoodXYZs().size() == MAXK);
      
    }else if(!ifs) {
      cerr << "file not exists. make fragmentgrid and save it " << fragment_grid_filename << endl;

      fragment_grids[i] = FragmentEnergyGrid(fragtemps[i], rotations_fragment, atom_grids, distance_grid, MAXK);
      assert(fragment_grids[i].getGoodXYZs().size() == MAXK);

      fragment_grids[i].writeFile(fragment_grid_filename);
    }else {
      cerr << "already exists. " << fragment_grid_filename << endl;

      fragment_grids[i] = FragmentEnergyGrid(i, fragment_grid_filename, MAXK);
      assert(fragment_grids[i].getGoodXYZs().size() == MAXK);
    }
    
    
    const vector<XYZ_param>& xyz_params = fragment_grids[i].getGoodXYZs();
    for(int j = 0; j < xyz_params.size(); j++) {
      const XYZ_param& xyz_param = xyz_params[j];
      Vector3d conv(fragment_grids[i].getGrid().convert(xyz_param.pos.x, xyz_param.pos.y, xyz_param.pos.z));
      cout << "f_smiles:" << fragsmiles[i] << ", f_xyz_param " << j << ":" << xyz_param.pos.x << "," << xyz_param.pos.y << "," << xyz_param.pos.z << "," << conv.x << "," << conv.y << "," << conv.z << "," << xyz_param.score << '\n';
    }
    cout << '\n';
    
  } 



  auto s2 = std::chrono::system_clock::now();
  fgrid_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  
  logs::lout << logs::info << "[TIME STAMP] fgrid_time : " << fgrid_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] END MAKE FRAGGRID" << endl;
  
  cout << "fgrid_time:" << fgrid_time.count() << endl;
  
  // フラグメントグリッドの最良スコアを集めた表の生成
  vector<pair<int, fltype>> fgrids_best_score = makeFgridsBestScore(fragment_grids);
  assert(fgrids_best_score.size() == ftmpsz);
  for(int i = 0; i < fgrids_best_score.size(); i++) {
    assert(fragment_grids[fgrids_best_score[i].first].getBestScore() == fgrids_best_score[i].second);
  }
  
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  cout << endl << endl;
  cout << "ligtypes,conftypes,fragtypes,fnums: " << lig_kind_sz << "," << ligs_sz << "," << ftmpsz << "," << fnums << endl;
  cout << endl << endl;
  if((config.per_frag == 10 and config.per_lig == 10)
  or (config.per_frag == 100 and config.per_lig == 100)) {
    printDatasetInfo(ligands_mol, fragtemps, fragvecs, same_ligs, frag_lig_list, fgrids_best_score, fragsmiles, lig_kind_sz, ligs_sz, ftmpsz, fnums); 
  }
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
//-----------------------
  
  int real_param_num=0;
  int tot_param_num=0;
  
  const int pf = config.per_frag;
  const int pl = config.per_lig;
  
  map<tuple<int, int, int>, int> mp;
  for(int isf1f2 = 0; isf1f2 <= 1; isf1f2++) {
  for(int pc = 100; pc <= 100; pc += 10) {
  
  if(isf1f2 == 0 and (pf != 10 or pl != 10 or pc != 100)) continue;
  
  //config.top_fgrid_num = k;
  config.f1_selection = isf1f2 ? true: false;
  config.f2_selection = isf1f2 ? true: false;
  config.per_frag = (fltype)pf / 100;
  config.per_lig  = (fltype)pl / 100;
  config.per_conf = (fltype)pc / 100;
  
  
  
  logs::lout << endl << endl;
  cout << endl << endl;
  
  std::chrono::milliseconds makeF1_time(0);
  std::chrono::milliseconds frag_selection_time(0);
  std::chrono::milliseconds lig_selection_time(0);
  std::chrono::milliseconds makeF2_time(0);
  std::chrono::milliseconds conf_selection_time(0);
  std::chrono::milliseconds detailed_score_calc_time(0);
  std::chrono::milliseconds real_time(0);
  
  auto r1 = std::chrono::system_clock::now();
  
  logs::lout << logs::info << "[TIME STAMP] START FRAG SELECTION" << endl;
  
  // 2. フラグメントの選抜
  vector<int> original_frags(ftmpsz);
  for(int i = 0; i < original_frags.size(); i++) {
    original_frags[i] = fgrids_best_score[i].first;
  }
  
  vector<int> selected_frags;
  if(config.f1_selection) {
    s1 = std::chrono::system_clock::now();
  
    selected_frags = frag_selection(fgrids_best_score, utils::round(ftmpsz * config.per_frag), config.thre_frag);

    s2 = std::chrono::system_clock::now();
    frag_selection_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  }else {
    selected_frags = original_frags;
  }
  
  
  logs::lout << logs::info << "[TIME STAMP] frag_selection_time : " << frag_selection_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] END FRAG SELECTION" << endl;
//-----------------------
  logs::lout << logs::info << "[TIME STAMP] START F1 SELECTION" << endl;
  
  // 3. F1選抜

  // 3.1 F1の生成
  s1 = std::chrono::system_clock::now();
  
  vector<pair<int, fltype>> f1 = makeF1(fgrids_best_score, selected_frags.size()); 

  s2 = std::chrono::system_clock::now();
  makeF1_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  
  for(int i = 0; i < f1.size(); i++) {
    assert(fragment_grids[f1[i].first].getBestScore() == f1[i].second);
  }


  vector<int> selected_ligs;
  vector<pair<int, fltype>> selected_ligs_score({});
  if (config.f1_selection) {
    // 3.2 化合物の選抜
    s1 = std::chrono::system_clock::now();

    selected_ligs_score = lig_selection(f1, frag_lig_list, original_ligs, utils::round(original_ligs.size() * config.per_lig), config.thre_lig);
    
    for(int i = 0; i < selected_ligs_score.size(); i++){
      selected_ligs.push_back(selected_ligs_score[i].first);
    }
    
    s2 = std::chrono::system_clock::now();
    lig_selection_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  }else {
    selected_ligs = original_ligs;
  }
  
  logs::lout << logs::info << "[TIME STAMP] makeF1_time : " << makeF1_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] lig_selection_time : " << lig_selection_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] END F1 SELECTION" << endl;
//-----------------------
  logs::lout << logs::info << "[TIME STAMP] START F2 SELECTION" << endl;
  
  // 4. F2選抜

  // 化合物レベル -> 配座レベルに落とす
  vector<int> original_confs;
  for(int i = 0; i < selected_ligs.size(); i++) {
    int lig_ind = selected_ligs[i];

    for(auto conf_ind: same_ligs[lig_ind]) {
      original_confs.push_back(conf_ind);
    }
  }
  
  
  vector<bool> is_selefrag(ftmpsz, false);
  for(int i = 0; i < selected_frags.size(); i++) {
    is_selefrag[selected_frags[i]] = true;  
  }
  
  set<pair<int, int>> f2pair;
  for(int i = 0; i < original_confs.size(); i++) {
    int conf_ind = original_confs[i];

    FragmentsVector& fv = fragvecs[conf_ind];
    int sz = fv.size();
    for(int j = 0; j < sz; j++) {
      int fj_id = fv.getvec(j).fragid;
      
      if(!is_selefrag[fj_id]) continue;
      
      for(int k = j + 1; k < sz; k++) {
        int fk_id = fv.getvec(k).fragid;
        
        if(!is_selefrag[fk_id]) continue;
        
        int id0 = std::min(fj_id, fk_id);
        int id1 = std::max(fj_id, fk_id);
          
        f2pair.insert(make_pair(id0, id1));  
        
      }
    }
  }
  //cout << "f2pair: " << f2pair.size() << endl;
  
  vector<int> selected_confs;
  vector<pair<int, fltype>> selected_confs_score({});
  if (config.f2_selection) {
    // 4.1 F2の生成
    s1 = std::chrono::system_clock::now();
    
    const int DMAX = utils::round(Vector3d(score_num.x, score_num.y, score_num.z).abs() + 1);
  
    map<pair<int, int>, vector<fltype>> f2 = makeF2(fragment_grids, f2pair, config.top_fgrid_num, DMAX); // 立方体の対角線が取り得る長さの最大
    
    s2 = std::chrono::system_clock::now();
    makeF2_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);

    // 4.2 配座の選抜
    s1 = std::chrono::system_clock::now();
    
    selected_confs_score = conf_selection(f1, f2, fragvecs, DMAX, config.grid.score_pitch_x, original_confs, std::min((int)original_confs.size(), utils::round(ligs_sz * config.per_conf)), config.thre_conf);
    
    for(int i = 0; i < selected_confs_score.size(); i++){
      selected_confs.push_back(selected_confs_score[i].first);
    }
    
    s2 = std::chrono::system_clock::now();
    conf_selection_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  }else {
    selected_confs = original_confs;
  }

  logs::lout << logs::info << "[TIME STAMP] makeF2_time : " << makeF2_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] conf_selection_time : " << conf_selection_time.count() << endl;
  logs::lout << logs::info << "[TIME STAMP] END F2 SELECTION" << endl;
//-----------------------
  logs::lout << logs::info << "[TIME STAMP] START CALCULATING BY FRAGGRID" << endl;
  
  
  if(config.f1_selection == false and config.f2_selection == false) {
    logs::lout << "not select at all" << endl;
    cout << "not select at all" << endl;
  }else {
    logs::lout << tot_param_num << " th param set" << endl;    
    cout << tot_param_num << " th param set" << endl;
    tot_param_num++;
  }

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  vector<vector<int>> original_same_ligs(lig_kind_sz, vector<int>());
  for(int i = 0; i < original_confs.size(); i++) {
    int conf_ind = original_confs[i];
    const string& identifier = ligands_mol[conf_ind].getIdentifier();
    original_same_ligs[lig_map[identifier]].push_back(conf_ind);
  }
  
  vector<vector<int>> selected_same_ligs(lig_kind_sz, vector<int>());
  for(int i = 0; i < selected_confs.size(); i++) {
    int conf_ind = selected_confs[i];
    const string& identifier = ligands_mol[conf_ind].getIdentifier();
    selected_same_ligs[lig_map[identifier]].push_back(conf_ind);
  }
  
  for(int i = 0; i < lig_kind_sz; i++) {
    cout << "conf transition" << setw(10) << i << ": "
    << setw(10) << same_ligs[i].size() << ", "
    << setw(10) << original_same_ligs[i].size() << ", "
    << setw(10) << selected_same_ligs[i].size() << "\n";
  }
  
  if(isf1f2 == 0) {
    cout << "no params" << endl;
  }else {
    cout << "params (pf,pl,pc, k)                     :" << pf << "," << pl << "," << pc << ", " << config.top_fgrid_num << endl;
  }


  cout << "origin,percent,threshold (lig_selection) :" << original_ligs.size() << "," <<  utils::round(original_ligs.size() * config.per_lig) << "," << selected_ligs.size() << endl;


  cout << "origin,percent,threshold (conf_selection):" << original_confs.size() << "," << std::min((int)original_confs.size(), utils::round(ligs_sz * config.per_conf)) << "," << selected_confs.size() << endl;  
  
  
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  if(isf1f2 == 0) {
    ;  
  }else {
    
    tuple<int, int, int> kumi = make_tuple(pf, pl, selected_confs.size());
    if(mp.count(kumi)) {
      continue;
    }else {
      mp[kumi] = 1;
      real_param_num++;
    } 
    
  }
  
  
  // 5. 詳細なスコア計算
  vector<pair<int, fltype>> fgrid_ranking(lig_kind_sz);
  for(int i = 0; i < lig_kind_sz; i++) {
    fgrid_ranking[i] = make_pair(i, LIMIT_ENERGY);
  }
  
  s1 = std::chrono::system_clock::now();
  
  for (int i = 0; i < selected_confs.size(); ++i) { //conformer loop
    int lig_ind = selected_confs[i];
    const Molecule& mol = ligands_mol[lig_ind];
    const string& identifier = mol.getIdentifier();
    // logs::lout << logs::info << (lig_ind + 1) << "th ligand : " << title << endl;

    int rotsz = rotations_ligand.size();
    int sz = fragvecs[lig_ind].size(); //化合物が持つフラグメントの数
    //各回転方向における化合物の重心に対するフラグメントの重心の相対位置
    vector<vector<int> > dx(rotsz, vector<int>(sz));
    vector<vector<int> > dy(rotsz, vector<int>(sz));
    vector<vector<int> > dz(rotsz, vector<int>(sz));

    vector<EnergyGrid> scores(rotsz, EnergyGrid(atom_grids[0].getCenter(), search_pitch, search_num, mol.getIntraEnergy()));
    // vector<EnergyGrid> scores(rotsz, EnergyGrid(atom_grids[0].getCenter(), search_pitch, search_num, 0.0));

    for (int rotid = 0; rotid < rotsz; ++rotid) {
      FragmentsVector fv = fragvecs[lig_ind];
      fv.rotate(rotations_ligand, rotid);
      //化合物の重心に対するフラグメントの重心の相対位置を求め保存
      for (int j = 0; j < sz; ++j) {
        dx[rotid][j] = utils::round(fv.getvec(j).pos.x / config.grid.score_pitch_x); 
        dy[rotid][j] = utils::round(fv.getvec(j).pos.y / config.grid.score_pitch_y);
        dz[rotid][j] = utils::round(fv.getvec(j).pos.z / config.grid.score_pitch_z);
      }
    }

    //auto t1 = std::chrono::system_clock::now();

    for (int j = 0; j < sz; ++j) { //fragment loop
      int fragid = fragvecs[lig_ind].getvec(j).fragid;
      const FragmentEnergyGrid& fg = fragment_grids[fragid];

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
    //auto t2 = std::chrono::system_clock::now();

    //fgrid_time += std::chrono::duration_cast< std::chrono::milliseconds >(t2 - t1);

    fltype best_conf_score = LIMIT_ENERGY;
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
              if (best_conf_score > score) {
                best_conf_score = score;
              }
            }
          }
        }
      }
    }
    
    fltype best_lig_score = fgrid_ranking[ind].second;
    if (best_lig_score > best_conf_score) {
      fgrid_ranking[ind].second = best_conf_score;
    }
    //auto t3 = std::chrono::system_clock::now();
    //real_time += std::chrono::duration_cast< std::chrono::milliseconds >(t3 - t2);

    // logs::lout << logs::info << "end calc real score" << endl;
  }
  
  s2 = std::chrono::system_clock::now();
  detailed_score_calc_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  

  logs::lout << logs::info << "[TIME STAMP] END CALCULATING BY FRAGGRID" << endl;
  logs::lout << logs::info << "[TIME STAMP] detailed_score_calc_time : " << detailed_score_calc_time.count() << endl;
  //logs::lout << logs::info << "[TIME STAMP] realcalc_time : " << real_time.count() << endl;

  logs::lout << logs::info << "[TIME STAMP] START RANKING BY FRAGGRID" << endl;

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  
  map<int, bool> ex;
  vector<pair<int, fltype>> lig_level_sele_confs_score;
  
  for(int i = 0; i < selected_confs_score.size(); i++) {
    int conf_ind = selected_confs_score[i].first;
    const Molecule& mol = ligands_mol[conf_ind];
    const string& identifier = mol.getIdentifier();
    int ind = lig_map[identifier];
    
    if(ex.count(ind)) {
      ;
    }else {
      ex[ind] = true;
      lig_level_sele_confs_score.push_back(make_pair(ind,selected_confs_score[i].second));
    }
    
  }
  
  assert(fgrid_ranking.size() == lig_kind_sz);
  sort(fgrid_ranking.begin(), fgrid_ranking.end(), [&](const pair<int, fltype>& p, const pair<int, fltype>& q){ return p.second < q.second; });
  
  //---
  
  cout << "sele lig ranking:"; // 化合物選抜後の化合物ランキング
  for(int i = 0; i < selected_ligs_score.size(); i++) {
    cout <<  selected_ligs_score[i].first << ",";
  }
  cout << endl;
  
  cout << "sele lig score:"; // 化合物選抜後の化合物スコア
  for(int i = 0; i < selected_ligs_score.size(); i++) {
    cout <<  selected_ligs_score[i].second << ",";
  }
  cout << endl;
  
  //---
  
  cout << "sele conf on conf level ranking:"; // 配座選抜後の配座ランキング(配座レベルで化合物IDを出力)
  for(int i = 0; i < selected_confs_score.size(); i++) {
    int conf_ind = selected_confs_score[i].first;
    const Molecule& mol = ligands_mol[conf_ind];
    const string& identifier = mol.getIdentifier();
    int ind = lig_map[identifier];
    cout <<  ind << ",";
  }
  cout << endl;
  
  cout << "sele conf on conf level score:"; // 配座選抜後の配座スコア
  for(int i = 0; i < selected_confs_score.size(); i++) {
    cout <<  selected_confs_score[i].second << ",";
  }
  cout << endl;
  
  //---
  
  cout << "sele conf on lig level ranking:"; // 配座選抜後の"化合物"ランキング
  for(int i = 0; i < lig_level_sele_confs_score.size(); i++) {
    cout <<  lig_level_sele_confs_score[i].first << ",";
  }
  cout << endl;
  
  cout << "sele conf on lig level score:"; // 配座選抜後の"化合物"スコア
  for(int i = 0; i < lig_level_sele_confs_score.size(); i++) {
    cout <<  lig_level_sele_confs_score[i].second << ",";
  }
  cout << endl;
  
  //---

  cout << "detailed ranking:"; // detailed 後の化合物ランキング
  for (int i = 0; i < fgrid_ranking.size(); ++i) {
    cout << fgrid_ranking[i].first << ",";
  }
  cout << endl;
  
  cout << "detailed score:"; // detailed 後の化合物スコア
  for (int i = 0; i < fgrid_ranking.size(); ++i) { 
    cout << fgrid_ranking[i].second << ",";
  }
  cout << endl;
  
  //---
  
  cout << "time:";
  cout << frag_selection_time.count() + makeF1_time.count() + lig_selection_time.count() + makeF2_time.count() + conf_selection_time.count() + detailed_score_calc_time.count() << "=" << frag_selection_time.count() << "+" << makeF1_time.count() << "+" << lig_selection_time.count() << "+" << makeF2_time.count() << "+" << conf_selection_time.count()<< "+" << detailed_score_calc_time.count() << endl;
  
  auto r2 = std::chrono::system_clock::now();
  real_time = std::chrono::duration_cast< std::chrono::milliseconds >(r2 - r1);
  cout << "real time:" << real_time.count() << endl;




//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  
  
  logs::lout << logs::info << "[TIME STAMP] END RANKING BY FRAGGRID" << endl;
  
  if(!(isf1f2 == 0 and pf == 10 and pl == 10 and pc == 100)) continue;
  
///*  
  std::chrono::milliseconds local_opt_time(0);
  s1 = std::chrono::system_clock::now();
  
  //vector<pair<fltype, string>> ranking;
  vector<tuple<fltype, string, fltype>> ranking;
  ranking.reserve(lig_kind_sz);
  OpenBabel::outputOBMol outputs(config.output_file);
  OpenBabel::outputOBMol orig_outputs(config.output_file + "_before" + ".sdf");
  //ofstream outputcsv(config.output_file + "fraggrid__" + getDate() + ".csv");
  ofstream outputcsv(config.output_file + "_fraggrid" + ".csv");

  logs::lout << logs::info << "start pre-calculate energy" << endl;
  // for optimize
  // calc = EnergyCalculator(0.95, 0.0);

  // Optimizer opt(receptor_mol);
  // vector<AtomEnergyGrid> opt_atom_grids = AtomEnergyGrid::makeAtomGrids(atom_grids[0].getCenter(), atom_grids[0].getPitch(), atom_grids[0].getNum(), receptor_mol, calc);
  Optimizer_Grid opt_grid(atom_grids, gradient_grids);

  logs::lout << logs::info << "end pre-calculate energy" << endl;


  logs::lout << logs::info << "[TIME STAMP] START OPTIMIZING AND RANKING" << endl;

  int opt_cnt = 0;
  cout << endl << endl;
  Point3d<int> inner_width(config.grid.inner_width_x,
                           config.grid.inner_width_y,
                           config.grid.inner_width_z);
  
  
  for (const auto& p : lig_map) {
    const string& identifier = p.first;
    int i = p.second;
    const vector<pos_param>& poss = pos_param_vec[i].getValues();
    //vector<pair<fltype, pair<int, Molecule> > > out_mols(poss.size());
    cout << "opt count" << opt_cnt << "," << identifier << endl;
    opt_cnt++;
    cout<<"poss size " << poss.size() << endl;
    vector<pair<fltype, tuple<int, Molecule, fltype, Molecule> > > out_mols(poss.size());
    int good=0;
    for (int j = 0; j < (int)poss.size(); j++) {
      const pos_param& param = poss[j];
      Molecule mol = ligands_mol[param.inp_ind];
      mol.rotate(rotations_ligand[param.rotid]);
      Vector3d pos = search_grid.convert(param.x, param.y, param.z);
      mol.translate(pos);
      //fltype opt_score = opt_grid.optimize(mol);
      cout<<"title:"<<mol.gettitle()<<",pose:"<<j<<'\n';
      fltype opt_score = opt_grid.optimizeBySteepest(mol, config.max_iterations, config.step_size);

      Molecule orig_mol = ligands_mol[param.inp_ind];
      orig_mol.rotate(rotations_ligand[param.rotid]);
      orig_mol.translate(pos);
      fltype orig_score = 0.0;
      for (auto& a : orig_mol.getAtoms()) {
        orig_score += atom_grids[a.getXSType()].getEnergy(a); //Optimizer_Grid::calcscore(mol) と同値
      }
      good+=(orig_score>opt_score);
      // fltype opt_score = param.score;

      //out_mols[j] = make_pair(opt_score, make_pair(param.inp_ind, mol));
      out_mols[j] = make_pair(opt_score, make_tuple(param.inp_ind, mol, orig_score, orig_mol));
    }
    //cout<<"kaizen "<<good<<endl;
    sort(out_mols.begin(), out_mols.end()); // opt_score でソート
    fltype best_intra = LIMIT_ENERGY;
    fltype best_score = LIMIT_ENERGY;
    fltype best_orig_score = LIMIT_ENERGY;
    int nrots = -999;
    fltype row_score = LIMIT_ENERGY;
    fltype row_orig_score = LIMIT_ENERGY;
    if (out_mols.size() > 0) {
      //best_intra = ligands_mol[out_mols[0].second.first].getIntraEnergy();
      //best_score = out_mols[0].first / (1 + 0.05846 * ligands_mol[out_mols[0].second.first].getNrots());
      best_intra = ligands_mol[get<0>(out_mols[0].second)].getIntraEnergy();
      best_score = out_mols[0].first / (1 + 0.05846 * ligands_mol[get<0>(out_mols[0].second)].getNrots());
      best_orig_score = get<2>(out_mols[0].second) / (1 + 0.05846 * ligands_mol[get<0>(out_mols[0].second)].getNrots());
      nrots = ligands_mol[get<0>(out_mols[0].second)].getNrots();
      row_score = out_mols[0].first;
      row_orig_score = get<2>(out_mols[0].second);
    }
    //ranking.push_back(make_pair(best_score, identifier));
    ranking.push_back(make_tuple(best_score, identifier, best_orig_score));

    //outputcsv << identifier << "," << best_score << endl;
    outputcsv << identifier << "," << best_score << "," << best_orig_score << "," << nrots << "," << best_intra << "," << row_score << "," << row_orig_score << endl;

    // logs::lout << "mol : " << title << endl;

    for (int j = 0; j < top_num && j < out_mols.size(); ++j) {
      //int lig_ind = out_mols[j].second.first;
      int lig_ind = get<0>(out_mols[j].second);
      fltype score = (out_mols[j].first + ligands_mol[lig_ind].getIntraEnergy() - best_intra) / (1 + 0.05846 * ligands_mol[lig_ind].getNrots());
      // fltype score = (out_mols[j].first) / (1 + 0.05846 * ligands_mol[lig_ind].getNrots());
      logs::lout << "  " << (j + 1) << "th pose's score : " << score << endl;
      OpenBabel::OBMol obmol = ligands[lig_ind];
      //OpenBabel::UpdateCoords(obmol, out_mols[j].second.second);
      OpenBabel::UpdateCoords(obmol, get<1>(out_mols[j].second));
      outputs.write(obmol);

      OpenBabel::OBMol orig_obmol = ligands[lig_ind];
      OpenBabel::UpdateCoords(orig_obmol, get<3>(out_mols[j].second));
      orig_outputs.write(orig_obmol);
      // obmol.AddHydrogens();
      // Molecule mol = format::toFragmentMol(obmol);
    }

  }
  cout << endl << endl;
  outputs.close();
  outputcsv.close();
  orig_outputs.close();

  assert(ranking.size() == lig_kind_sz);
  sort(ranking.begin(), ranking.end());
  for (int i = 0; i < lig_kind_sz; ++i) {
    //logs::lout << (i + 1) << "th ligand : " << ranking[i].second << endl;
    logs::lout << (i + 1) << "th ligand : " << get<1>(ranking[i]) << endl;
    //logs::lout << "score : " << ranking[i].first << endl;
    logs::lout << "score : " << get<0>(ranking[i]) << endl;
  }

  logs::lout << logs::info << "[TIME STAMP] END OPTIMIZING AND RANKING" << endl;

  // logs::lout << "predict reduce cost : " << pred_reduce << endl;
  logs::lout << "real reduce cost    : " << reduces << endl;
  cerr << "real reduce cost    : " << reduces << endl;
  //assert(config.reuse_grid != format::DockingConfiguration::REUSE_OFFLINE || reduces == pred_reduce);

  // OpenBabel::outputOBMolsToSDF(config.output_file, outputobmols);

  logs::lout << logs::info << "################ Program end ################" << endl;

  logs::lout << logs::info << "[FINAL_RESULT] "
      << ligs_sz << "/" << config.getReuseGridString() << "/" << (config.reorder ? "reorder" : "no_reorder") << "/" << config.mem_size << "/" << config.grid.inner_width_x << ", "
      << "fragment types : " << ftmpsz << ", "
      << "fragment num : " << fnums << ", "
      << "all cost : " << allcost << ", "
      << "minimum cost : " << minimum_cost << ", "
      << "reduce cost : " << reduces << ", "
      << "real cost : " << allcost - reduces << ", "
      << "realcalc_time : " << real_time.count() << endl;

  // logs::lout << logs::info << "[RANKING_SCORE_RESULT]" << endl;
  // for (int i = 0; i < lig_kind_sz; ++i) {
  //   if (i > 0)
  //     logs::lout << " " << ranking[i].first;
  //   else
  //     logs::lout << ranking[i].first;
  // }
  // logs::lout << endl;
  
  cout << "REstretto ranking:";
  for (int i = 0; i < lig_kind_sz; ++i) {
    cout << lig_map[get<1>(ranking[i])] << ",";
  }
  cout << endl;
  
  cout << "REstretto score:";
  for (int i = 0; i < lig_kind_sz; ++i) {
    cout << get<0>(ranking[i]) << ",";
  }
  cout << endl;
  
  cout << "REstretto orig score:";
  for (int i = 0; i < lig_kind_sz; ++i) {
    cout << get<2>(ranking[i]) << ",";
  }
  cout << endl;
  
  s2 = std::chrono::system_clock::now();
  local_opt_time = std::chrono::duration_cast< std::chrono::milliseconds >(s2 - s1);
  cout << "local opt time:" << local_opt_time.count() << endl;
//*/  
  }}
  cout << endl;
  
  cout << "total num of param sets:" << tot_param_num <<endl;
  cout << "real num of param sets:" << real_param_num <<endl;
  cout << "edit check" << endl;
  cout << "##########" << endl << endl;
  logs::lout.close();
  return 0;
}
