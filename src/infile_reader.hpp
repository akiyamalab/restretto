#include "common.hpp"
#include "EnergyGrid.hpp"
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
    ReuseStrategy reuse_grid = ReuseStrategy::OFFLINE;
    bool reorder = true;
    int64_t mem_size;
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
