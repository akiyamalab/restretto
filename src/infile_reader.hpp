#include "common.hpp"
#include "InterEnergyGrid.hpp"
#include <string>
#include <vector>

#ifndef INFILE_READER_H_
#define INFILE_READER_H_

namespace format {
  struct SearchGrid {
    fragdock::Point3d<fltype> center, outer_width, inner_width, search_pitch, score_pitch;
  };

  struct DockingConfiguration {
    enum ReuseStrategy {OFFLINE, ONLINE, NONE};

    SearchGrid grid;
    std::vector<std::string> ligand_files;
    std::string receptor_file, output_file;
    std::string log_file, grid_folder, dxgrid_folder;
    std::string rotangs_file;
    ReuseStrategy reuse_grid = ReuseStrategy::OFFLINE;
    bool reorder = true;
    int64_t mem_size;
    int64_t poses_per_lig = 1;
    int64_t poses_per_lig_before_opt = 2000;
    fltype output_score_threshold = -3.0;
    fltype pose_min_rmsd = 0.5;
    bool no_local_opt = false;
    bool score_only = false;
    bool local_only = false;
    fltype local_max_rmsd = 1e10;
    fltype rad_scale = 0.95;
    const std::string getReuseGridString() {
      switch (reuse_grid) {
        case ReuseStrategy::OFFLINE: return "REUSE_OFFLINE";
        case ReuseStrategy::ONLINE : return "REUSE_ONLINE";
        case ReuseStrategy::NONE   : return "REUSE_NONE";
      }
      assert(0);
    }
    void checkConfigValidity() const;
  };

  DockingConfiguration ParseInFile(const char *filename);
} // namespace format

#endif
