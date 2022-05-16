#ifndef FRAGVEC_H_
#define FRAGVEC_H_

#include "common.hpp"
#include "Vector3d.hpp"

namespace fragdock {
  struct fragvec {
    Vector3d pos;
    // int rotid; //これをメンバに持っていないのはおそらく最良方向のスコアだけを用いる実装に現在なっているからであり，スコアの候補を増やす場合はリガンドの重心に対するフラグメントの重心(真上にある Vector3d pos)の情報に加えて，さらにその重心に対するフラグメントの回転方向 int rotid の情報が必要になる．
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

      std::sort(vecs.begin(), vecs.end()); //20行目にfragvecの大小が定義されているのでそれに従ってソートされる．
    }
    bool operator<(const FragmentsVector& o) const { return vecs < o.vecs; } //各要素を1つずつ大小比較していき差がついた時点でその差を大小の差と定義する．vector<pll>のソートでもやっていること．a=[2,4,4,11]とb=[2,4,5,9]なら3番目が初めて異なる要素であり4の方が5より小さいのでaが配列全体の大小としてbより小さいと定義される(初めて差がついた時点のみ考えるため，それ以降の要素の大小は関係ない)．vectorに関しての"<"演算子のオーバーロードがデフォルトでなされているからvecs < o.vecs と書くだけで良くなっている．．
  };
}
#endif
