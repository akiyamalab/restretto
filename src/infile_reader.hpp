#include "common.hpp"
#include <string>
#include <vector>

#ifndef INFILE_READER_H_
#define INFILE_READER_H_
namespace format {
  struct SearchGrid {
    fltype cx, cy, cz;
    fltype outer_width_x, outer_width_y, outer_width_z;
    fltype inner_width_x, inner_width_y, inner_width_z;
    fltype search_pitch_x, search_pitch_y, search_pitch_z;
    fltype score_pitch_x, score_pitch_y, score_pitch_z;
  };

  struct DockingConfiguration {
    static const int REUSE_OFFLINE = 1;
    static const int REUSE_ONLINE  = 0;
    static const int REUSE_NONE    = -1;

    SearchGrid grid;
    std::vector<std::string> ligand_files;
    std::string receptor_file, output_file;
    std::string log_file, grid_folder;
    std::string rotangs_file;
    int reuse_grid;
    bool reorder;
    int64_t mem_size;

    std::string fragment_grid_folder;
    bool f1_selection;
    bool f2_selection;
    int64_t top_fgrid_num;
    fltype per_frag;
    fltype per_lig;
    fltype per_conf;
    fltype thre_frag;
    fltype thre_lig;
    fltype thre_conf;

    int64_t max_iterations;
    fltype step_size;

    
    const std::string getReuseGridString() {
      switch (reuse_grid) {
        case REUSE_OFFLINE: return "REUSE_OFFLINE";
        case REUSE_ONLINE : return "REUSE_ONLINE";
        case REUSE_NONE   : return "REUSE_NONE";
      }
      assert(0);
    }
  };

  DockingConfiguration ParseInFile(const char *filename);
} // namespace format

#endif
