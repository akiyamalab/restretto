#include <vector>

#ifndef _MIN_VALUES_VECTOR_H_
#define _MIN_VALUES_VECTOR_H_

namespace utils{

  /**
   * @brief A class to store the minimum values of a vector.
  */
  template <typename T>
  class MinValuesVector {
    int size, threthold;
    std::vector<T> vec;

    void resize(int size) {
      std::nth_element(vec.begin(), vec.begin() + size, vec.end());
      vec.resize(size);
    }

  public:
    /* @param size The number how many elements will be stored. */
    MinValuesVector(int size) : size(size), threthold(size * 2) {}

    void push(const T& x) {
      vec.push_back(x);
      if (vec.size() >= threthold) {
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
