#ifndef FRAGVEC_H_
#define FRAGVEC_H_

#include "common.hpp"
#include "Vector3d.hpp"

namespace fragdock {
  /**
   * A reduced representation of a fragment.
  */
  struct fragvec {
    Vector3d pos; // center coordinate of a fragment
    int frag_idx; // index of a fragment
    int size;     // # of atoms in a fragment
    int rank;     // importance ranking of a fragment
    fragvec(const Vector3d& p, int fragid, int size) : frag_idx(fragid), size(size) {
      pos.x = p.x;
      pos.y = p.y;
      pos.z = p.z;
    }

    // set importance ranking of a fragment
    void setrank(int r) { rank = r; }
    bool operator<(const fragvec& o) const { return rank < o.rank; }
  };

  /**
   * Fragments which are decomposed from a molecule. 
  */
  class FragmentsVector {
    std::vector<fragvec> vecs;
  public:
    // # of fragments
    int size() const { return vecs.size(); }
    // Add a fragment
    void append(const fragvec& v) { vecs.push_back(v); }
    const std::vector<fragvec>& getvecs() const { return vecs; }
    const fragvec& getvec(int i) const { return vecs[i]; }
    void rotate(const std::vector<Vector3d>& rotations, int rotid) {
      for (auto& fv : vecs) {
        fv.pos.rotate(rotations[rotid]);
        // fv.rotid = RotMatrix[fv.rotid][rotid];
      }
    }

    // Translate all fragments
    void translate(const Vector3d& v) {
      for (auto& fv : vecs)
        fv.pos += v;
    }

    /** 
     * Sort fragments by (importance) rank
     * fragrank: fragment idx -> rank
    */
    void sort(const std::vector<int>& fragrank) {
      for (auto& fv : vecs)
        fv.setrank(fragrank[fv.frag_idx]);
      std::sort(vecs.begin(), vecs.end());
    }
    bool operator<(const FragmentsVector& o) const { return vecs < o.vecs; }
  };
}
#endif
