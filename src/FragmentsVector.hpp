#ifndef FRAGVEC_H_
#define FRAGVEC_H_

#include "common.hpp"
#include "Vector3d.hpp"

namespace fragdock {
  struct fragvec {
    Vector3d pos;
    // int rotid;
    int fragid; // equal temp_id
    int size;
    int rank;
    fragvec(const Vector3d& p, int fragid, int size) : fragid(fragid), size(size) {
      pos.x = p.x;
      pos.y = p.y;
      pos.z = p.z;
    }
    void setrank(int r) { rank = r; }
    bool operator<(const fragvec& o) const { return rank < o.rank; }
  };
  class FragmentsVector {
    std::vector<fragvec> vecs;
  public:
    int size() const { return vecs.size(); }
    void append(const fragvec& v) { vecs.push_back(v); }
    const std::vector<fragvec>& getvecs() const { return vecs; }
    const fragvec& getvec(int i) const { return vecs[i]; }
    void rotate(const std::vector<Vector3d>& rotations, int rotid) {
      for (auto& fv : vecs) {
        fv.pos.rotate(rotations[rotid]);
        // fv.rotid = RotMatrix[fv.rotid][rotid];
      }
    }
    void translate(const Vector3d& v) {
      for (auto& fv : vecs)
        fv.pos += v;
    }
    void sort(const std::vector<int>& fragrank) {
      for (auto& fv : vecs)
        fv.setrank(fragrank[fv.fragid]);

      std::sort(vecs.begin(), vecs.end());
    }
    bool operator<(const FragmentsVector& o) const { return vecs < o.vecs; }
  };
}
#endif
