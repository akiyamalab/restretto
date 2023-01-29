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
    static const int REUSE_OFFLINE = 1;
    static const int REUSE_ONLINE  = 0;
    static const int REUSE_NONE    = -1;

    SearchGrid grid;
    std::vector<std::string> ligand_files;
    std::string receptor_file, output_file;
    std::string log_file, grid_folder;
    std::string rotangs_file;
    int reuse_grid = REUSE_OFFLINE;
    bool reorder = true;
    int64_t mem_size;
    const std::string getReuseGridString() {
      switch (reuse_grid) {
        case REUSE_OFFLINE: return "REUSE_OFFLINE";
        case REUSE_ONLINE : return "REUSE_ONLINE";
        case REUSE_NONE   : return "REUSE_NONE";
      }
      assert(0);
    }
    void checkConfigValidity() const;
  };

  DockingConfiguration ParseInFile(const char *filename);
} // namespace format

#endif
