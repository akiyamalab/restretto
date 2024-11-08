#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include "infile_reader.hpp"


namespace format {
  DockingConfiguration ParseInFile(const char *filename){
    std::ifstream ifs(filename);
    if (ifs.fail()){
      logs::lerr << logs::error << "opening config file failed: " << filename << std::endl;
      abort();
    }
    DockingConfiguration conf;
    std::string buffer;
    while(!ifs.eof()){
      std::getline(ifs, buffer);
      if (boost::algorithm::starts_with(buffer, "OUTERBOX ")) {
        std::string str = buffer.substr(9);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.outer_width = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "INNERBOX ")) {
        std::string str = buffer.substr(9);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.inner_width = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "BOX_CENTER ")) {
        std::string str = buffer.substr(11);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.center = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "SEARCH_PITCH ")) {
        std::string str = buffer.substr(13);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.search_pitch = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "SCORING_PITCH ")) {
        std::string str = buffer.substr(14);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.score_pitch = fragdock::Point3d<fltype>(
          boost::lexical_cast<fltype>(vals[0]),
          boost::lexical_cast<fltype>(vals[1]),
          boost::lexical_cast<fltype>(vals[2]));
      }
      else if (boost::algorithm::starts_with(buffer, "REUSE_FRAG_GRID ")) {
        std::string str = buffer.substr(16);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str == "false" || str == "none" || str == "noreuse") {
          conf.reuse_grid = DockingConfiguration::ReuseStrategy::NONE;
        }
        else if (str == "online") {
          conf.reuse_grid = DockingConfiguration::ReuseStrategy::ONLINE;
        }
        else {
          conf.reuse_grid = DockingConfiguration::ReuseStrategy::OFFLINE;
        }
      }
      else if (boost::algorithm::starts_with(buffer, "REORDER_LIGANDS ")) {
        std::string str = buffer.substr(16);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str == "true") conf.reorder = true;
        else if (str == "false") conf.reorder = false;
        else {
          logs::lerr << logs::error << "invalid value for REORDER_LIGANDS. It must be either 'true' or 'false'." << std::endl;
          abort();
        }
      }
      else if (boost::algorithm::starts_with(buffer, "MEMORY_SIZE ")) {
        // fraggrid_main の方でデフォルト値を設定しちゃったから使えないような気がする
        std::string str = buffer.substr(12);
        boost::algorithm::trim(str);
        conf.mem_size = boost::lexical_cast<int64_t>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "RECEPTOR ")) {
        conf.receptor_file = buffer.substr(9);
      }
      else if (boost::algorithm::starts_with(buffer, "LIGAND ")) {
        conf.ligand_files.push_back(buffer.substr(7));
      }
      else if (boost::algorithm::starts_with(buffer, "OUTPUT ")) {
        conf.output_file = buffer.substr(7);
      }
      // if(boost::algorithm::starts_with(buffer, "FRAGMENT ")){
      //   conf.fragment_file = buffer.substr(9);
      // }
      else if (boost::algorithm::starts_with(buffer, "LOG ")) {
        conf.log_file = buffer.substr(4);
      }
      // else if (boost::algorithm::starts_with(buffer, "FORCEFIELD ")) {
      //   conf.forcefield_file = buffer.substr(11);
      // }
      else if (boost::algorithm::starts_with(buffer, "GRID_FOLDER ")) {
        conf.grid_folder = buffer.substr(12);
      }
      else if (boost::algorithm::starts_with(buffer, "ROTANGS ")) {
        conf.rotangs_file = buffer.substr(8);
      }
      else if (boost::algorithm::starts_with(buffer, "POSES_PER_LIG ")) {
        std::string str = buffer.substr(14);
        boost::algorithm::trim(str);
        conf.poses_per_lig = boost::lexical_cast<int64_t>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "POSES_PER_LIG_BEFORE_OPT ")) {
        std::string str = buffer.substr(25);
        boost::algorithm::trim(str);
        conf.poses_per_lig_before_opt = boost::lexical_cast<int64_t>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "OUTPUT_SCORE_THRESHOLD ")) {
        std::string str = buffer.substr(23);
        boost::algorithm::trim(str);
        conf.output_score_threshold = boost::lexical_cast<fltype>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "MIN_RMSD ")) {
        std::string str = buffer.substr(9);
        boost::algorithm::trim(str);
        conf.pose_min_rmsd = boost::lexical_cast<fltype>(str);
      }
      else if (boost::algorithm::starts_with(buffer, "NO_LOCAL_OPT ")) {
        std::string str = buffer.substr(13);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str == "true") conf.no_local_opt = true;
        else if (str == "false") conf.no_local_opt = false;
        else {
          logs::lerr << logs::error << "Error: invalid value for NO_LOCAL_OPT. It must be either 'true' or 'false'." << std::endl;
          abort();
        }
      }
      else if (boost::algorithm::starts_with(buffer, "SCORE_ONLY ")) {
        std::string str = buffer.substr(11);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str == "true") conf.score_only = true;
        else if (str == "false") conf.score_only = false;
        else {
          logs::lerr << logs::error << "Error: invalid value for SCORE_ONLY. It must be either 'true' or 'false'." << std::endl;
          abort();
        }
      }
      else if (boost::algorithm::starts_with(buffer, "LOCAL_ONLY ")) {
        std::string str = buffer.substr(11);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str == "true") conf.local_only = true;
        else if (str == "false") conf.local_only = false;
        else {
          logs::lerr << logs::error << "Error: invalid value for LOCAL_ONLY. It must be either 'true' or 'false'." << std::endl;
          abort();
        }
      }
      else if (boost::algorithm::starts_with(buffer, "LOCAL_MAX_RMSD ")) {
        std::string str = buffer.substr(15);
        boost::algorithm::trim(str);
        conf.local_max_rmsd = boost::lexical_cast<fltype>(str);
      }
    }

    return conf;
  }
  void DockingConfiguration::checkConfigValidity() const {
    // ratio of search_pitch and score_pitch must be integer
    const fragdock::Point3d<int> ratio = utils::round(grid.search_pitch / grid.score_pitch);
    assert(abs(grid.score_pitch.x * ratio.x - grid.search_pitch.x) < EPS);
    assert(abs(grid.score_pitch.y * ratio.y - grid.search_pitch.y) < EPS);
    assert(abs(grid.score_pitch.z * ratio.z - grid.search_pitch.z) < EPS);

    // memory size must be enough to store at least ONE grid
    const fragdock::Point3d<int> num = utils::ceili(grid.outer_width / 2 / grid.score_pitch) * 2 + 1; // # of grid points per axis
    int FGRID_SIZE = (int)((mem_size * 1024 * 1024) / ((int64_t) num.x * num.y * num.z * sizeof(fltype))); // # of grids that can be stored in memory
    assert(FGRID_SIZE > 0);

    // option conflict check (do not set at the same time)
    if (score_only && local_only) {
      logs::lerr << logs::error << "'score-only' and 'local-only' cannot be true at the same time." << std::endl;
      abort();
    }
    if (score_only && no_local_opt) {
      logs::lerr << logs::warn << "'score-only' and 'no-local-opt' are both true. This means the same as 'score-only' only." << std::endl;
    }
    if (local_only && no_local_opt) {
      logs::lerr << logs::warn << "Warning: 'local-only' and 'no-local-opt' are both true. This may the same as 'score-only'." << std::endl;
    }

  }
} // namespace format
