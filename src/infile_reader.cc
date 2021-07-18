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
      std::cerr << "opening grid file failed:" << filename << std::endl;
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
        conf.grid.outer_width_x = boost::lexical_cast<fltype>(vals[0]);
        conf.grid.outer_width_y = boost::lexical_cast<fltype>(vals[1]);
        conf.grid.outer_width_z = boost::lexical_cast<fltype>(vals[2]);
      }
      else if (boost::algorithm::starts_with(buffer, "INNERBOX ")) {
        std::string str = buffer.substr(9);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.inner_width_x = boost::lexical_cast<fltype>(vals[0]);
        conf.grid.inner_width_y = boost::lexical_cast<fltype>(vals[1]);
        conf.grid.inner_width_z = boost::lexical_cast<fltype>(vals[2]);
      }
      else if (boost::algorithm::starts_with(buffer, "BOX_CENTER ")) {
        std::string str = buffer.substr(11);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.cx = boost::lexical_cast<fltype>(vals[0]);
        conf.grid.cy = boost::lexical_cast<fltype>(vals[1]);
        conf.grid.cz = boost::lexical_cast<fltype>(vals[2]);
      }
      else if (boost::algorithm::starts_with(buffer, "SEARCH_PITCH ")) {
        std::string str = buffer.substr(13);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.search_pitch_x = boost::lexical_cast<fltype>(vals[0]);
        conf.grid.search_pitch_y = boost::lexical_cast<fltype>(vals[1]);
        conf.grid.search_pitch_z = boost::lexical_cast<fltype>(vals[2]);
      }
      else if (boost::algorithm::starts_with(buffer, "SCORING_PITCH ")) {
        std::string str = buffer.substr(14);
        std::vector<std::string> vals;
        boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
        boost::algorithm::trim(vals[0]);
        boost::algorithm::trim(vals[1]);
        boost::algorithm::trim(vals[2]);
        conf.grid.score_pitch_x = boost::lexical_cast<fltype>(vals[0]);
        conf.grid.score_pitch_y = boost::lexical_cast<fltype>(vals[1]);
        conf.grid.score_pitch_z = boost::lexical_cast<fltype>(vals[2]);
      }
      else if (boost::algorithm::starts_with(buffer, "REUSE_FRAG_GRID ")) {
        std::string str = buffer.substr(16);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        if (str == "false" || str == "none" || str == "noreuse") {
          conf.reuse_grid = DockingConfiguration::REUSE_NONE;
        }
        else if (str == "online") {
          conf.reuse_grid = DockingConfiguration::REUSE_ONLINE;
        }
        else {
          conf.reuse_grid = DockingConfiguration::REUSE_OFFLINE;
        }
      }
      else if (boost::algorithm::starts_with(buffer, "REORDER_LIGANDS ")) {
        std::string str = buffer.substr(16);
        boost::algorithm::trim(str);
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        conf.reorder = (str != "false");
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
    }
    return conf;
  }
} // namespace format
