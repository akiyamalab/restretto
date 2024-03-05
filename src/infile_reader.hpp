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
    std::string log_file, grid_folder;
    std::string rotangs_file;
    std::string fragment_file;
    ReuseStrategy reuse_grid = ReuseStrategy::OFFLINE;
    bool reorder = true;
    int64_t mem_size;
    int64_t poses_per_lig = 1;
    fltype pose_rmsd = 0.5;
    bool no_local_opt = false;
    int capping_atomic_num, max_ring_size;
    bool do_carbon_capping, insert_fragment_id_to_isotope, merge_solitary;
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
