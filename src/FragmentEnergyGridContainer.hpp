#include "common.hpp"
#include "infile_reader.hpp"
#include "FragmentEnergyGrid.hpp"


#ifndef _FRAGMENT_ENERGY_GRID_CONTAINER_H_
#define _FRAGMENT_ENERGY_GRID_CONTAINER_H_

typedef format::DockingConfiguration::ReuseStrategy Strategy;

namespace fragdock {
  class FragmentEnergyGridContainer {
  private:
    int size;
    std::vector<FragmentEnergyGrid> grids;
    format::DockingConfiguration::ReuseStrategy strategy;
    std::vector<int> indices_to_save; // only for OFFLINE strategy
    int step = 0;
    std::vector<int> last_used;       // only for ONLINE strategy
    void checkValidity() { // only for OFFLINE strategy
      for (auto i : indices_to_save) {
        if (i >= size) {
          logs::lout << logs::error << "Invalid index to save: " << i << std::endl;
          throw std::runtime_error("Invalid index to save");
        }
      }
    }
    int search(const int fragid) const { // sequential search
      for (int i = 0; i < size; i++) {
        if (grids[i].frag_idx == fragid) return i;
      }
      return -1;
    }
    int find_lru_idx() const { // find least recently used
      int last = 1e9;
      int lru_idx = -1;
      for (int k = 0; k < size; ++k) {
        if (last_used[k] < last) {
          last = last_used[k];
          lru_idx = k;
        }
      }
      return lru_idx;
    }
  public:
    FragmentEnergyGridContainer() : size(0) {}
    FragmentEnergyGridContainer(const int size) : size(size) {
      grids = std::vector<FragmentEnergyGrid>(size);
      strategy = Strategy::ONLINE;
      last_used = std::vector<int>(size, -1);
    }
    FragmentEnergyGridContainer(const int size, const std::vector<int> &indices_to_save) : size(size), indices_to_save(indices_to_save) {
      grids = std::vector<FragmentEnergyGrid>(size);
      strategy = Strategy::OFFLINE;
      checkValidity();
    }
    FragmentEnergyGridContainer& operator=(const FragmentEnergyGridContainer& other) {
      if (this != &other) {
        size = other.size;
        grids = other.grids;
        strategy = other.strategy;
        indices_to_save = other.indices_to_save;
        step = other.step;
        last_used = other.last_used;
      }
      return *this;
    }
    ~FragmentEnergyGridContainer() {}
    void insert(const FragmentEnergyGrid& grid) {
      if (strategy == Strategy::ONLINE) {
        int lru_idx = find_lru_idx();
        grids[lru_idx] = grid;
        last_used[lru_idx] = step;
      } else { // strategy == Strategy::OFFLINE
        grids[indices_to_save[step]] = grid;
      }
    }
    bool isRegistered(const int fragid) const {
      if (strategy == Strategy::ONLINE) {
        return search(fragid) != -1;
      } else { // strategy == Strategy::OFFLINE
        return grids[indices_to_save[step]].frag_idx == fragid;
      } 
    }
    const FragmentEnergyGrid& get(const int fragid) {
      if (! isRegistered(fragid)) {
        logs::lout << logs::error << "FragmentEnergyGridContainer: FragmentEnergyGrid not found: " << fragid << std::endl;
        throw std::runtime_error("FragmentEnergyGridContainer: FragmentEnergyGrid not found");
      }

      if (strategy == Strategy::ONLINE) {
        int grid_idx = search(fragid);
        last_used[grid_idx] = step; // last used
        return grids[grid_idx];
      } else {
        return grids[indices_to_save[step]];
      }
    }
    void next() { step++; }
  };
}

#endif

    // for (int rotid = 0; rotid < rotsz; ++rotid) {
    //   FragmentsVector fv = fragvecs[lig_ind];
    //   fv.rotate(rotations_ligand, rotid);

    //   for (int j = 0; j < sz; ++j) {
    //     d[rotid][j] = round(fv.getvec(j).pos / config.grid.score_pitch);
    //   }
    // }

    // auto t1 = std::chrono::system_clock::now();

    // for (int j = 0; j < sz; ++j, ++frag_itr) {
    //   int fragid = fragvecs[lig_ind].getvec(j).frag_idx;
    //   int nextsp = -1;
    //   if (config.reuse_grid == format::DockingConfiguration::ReuseStrategy::OFFLINE || config.reuse_grid == format::DockingConfiguration::ReuseStrategy::NONE) {
    //     if (config.reuse_grid == format::DockingConfiguration::ReuseStrategy::OFFLINE)
    //       nextsp = nextgridsp[frag_itr];
    //     else
    //       nextsp = 0;

    //     if (fragment_grids[nextsp].frag_idx != fragid) {
    //       fragment_grids[nextsp] = FragmentEnergyGrid(frag_library[fragid], rotations_fragment, atom_grids, distance_grid);
    //     }
    //     else {
    //       reduces += fragvecs[lig_ind].getvec(j).size;
    //     }
    //   }
    //   else if (config.reuse_grid == format::DockingConfiguration::ReuseStrategy::ONLINE) {
    //     if (fgrid_ind[fragid] == -1) {
    //       int mi = frag_itr;
    //       for (int k = 0; k < FGRID_SIZE; ++k) {
    //         if (last_used[k] < mi) {
    //           mi = last_used[k];
    //           nextsp = k;
    //         }
    //       }
    //       assert(nextsp != -1);
    //       if (mi != -1) {
    //         // logs::lout << logs::info << "remove id : " << fragment_grids[nextsp].frag_idx << " fragment grid in vector at " << nextsp << endl;
    //         fgrid_ind[fragment_grids[nextsp].frag_idx] = -1;
    //       }
    //       // logs::lout << logs::info << "calc id : " << fragid << " fragment grid and store to vector at " << nextsp << endl;
    //       fragment_grids[nextsp] = FragmentEnergyGrid(frag_library[fragid], rotations_fragment, atom_grids, distance_grid);
    //       fgrid_ind[fragid] = nextsp;
    //     }
    //     else {
    //       nextsp = fgrid_ind[fragid];
    //       // logs::lout << logs::info << "get id : " << fragid << " fragment grid from vector at " << nextsp << endl;
    //       reduces += fragvecs[lig_ind].getvec(j).size;
    //     }
    //     last_used[fgrid_ind[fragid]] = frag_itr;
    //   }
    //   else {
    //     assert(0);
    //   }

    //   assert(fragment_grids[nextsp].frag_idx == fragid);
    //   const FragmentEnergyGrid& fg = fragment_grids[nextsp];

    //   #pragma omp parallel for // Calculation among rotation is independent
    //   for (int rotid = 0; rotid < rotsz; ++rotid) {
    //     // int rid = RotMatrix[fragvecs[lig_ind].getvec(j).rotid][rotid];
    //     Point3d<int> gs = to_score_num(0, score_num, search_num, ratio);
    //     for (int x = 0, gx = gs.x; x < search_num.x; ++x, gx += ratio.x)
    //     for (int y = 0, gy = gs.y; y < search_num.y; ++y, gy += ratio.y)
    //     for (int z = 0, gz = gs.z; z < search_num.z; ++z, gz += ratio.z)
    //       scores[rotid].addEnergy(x, y, z, fg.getGrid().getEnergy(gx + d[rotid][j].x, gy + d[rotid][j].y, gz + d[rotid][j].z));