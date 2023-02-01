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
    Strategy strategy;
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
