#include <vector>

#ifndef _MIN_VALUES_VECTOR_H_
#define _MIN_VALUES_VECTOR_H_

namespace utils{

  template <typename T>
  class MinValuesVector {
    int size, threthold;
    std::vector<T> vec;

    void resize(int size) { // ユーザは直接使わない．
      std::nth_element(vec.begin(), vec.begin() + size, vec.end()); //vec.begin()+size番目に小さい値よりも小さい値がvec.begin()+size番目に小さい値の前にくるように並び替える．
      vec.resize(size); // MinValuesVectorの resize ではなくstd::vector の resize なので注意．一つ上の処理後のこの処理により，Top-sizeのスコアはそのままに保ったまま vector の大きさをを size に縮小することが出来た．
    }

  public:
    MinValuesVector(int size) : size(size), threthold(size * 2) {}

    void push(const T& x) {
      vec.push_back(x);
      if (vec.size() >= threthold) { //vec.size()がsize+1となる度に逐一resizeをするのは時間の無駄？だから threthold = size * 2 というように多少余裕を持たせている？
        resize(size);
      }
    }

    const std::vector<T>& getValues() {
      if (vec.size() > size) {
        resize(size);
      }
      return vec;
    }
  };

}

#endif
