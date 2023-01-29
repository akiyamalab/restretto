#include "common.hpp"
#include "utils.hpp"

#ifndef _POINT3D_H_
#define _POINT3D_H_

namespace fragdock {
  template<class T>
  struct Point3d {
    T x, y, z;
    Point3d() : x(T()), y(T()), z(T()) {}
    Point3d(const T &x, const T &y, const T &z) : x(x), y(y), z(z) {}

    template<class U>
    auto operator+(const Point3d<U>& o) const -> Point3d<decltype(x+o.x)> {
      return Point3d<decltype(x+o.x)>(x+o.x, y+o.y, z+o.z);
    }
    auto operator+(const T& o) const -> Point3d<decltype(x+o)> {
      return Point3d<decltype(x+o)>(x+o, y+o, z+o);
    }
    template<class U>
    auto operator-(const Point3d<U>& o) const -> Point3d<decltype(x-o.x)> {
      return Point3d<decltype(x-o.x)>(x-o.x, y-o.y, z-o.z);
    }
    auto operator-(const T& o) const -> Point3d<decltype(x-o)> {
      return Point3d<decltype(x-o)>(x-o, y-o, z-o);
    }
    template<class U>
    auto operator*(const Point3d<U>& o) const -> Point3d<decltype(x*o.x)> {
      return Point3d<decltype(x*o.x)>(x*o.x, y*o.y, z*o.z);
    }
    auto operator*(const T& o) const -> Point3d<decltype(x*o)> {
      return Point3d<decltype(x*o)>(x*o, y*o, z*o);
    }
    template<class U>
    auto operator/(const Point3d<U>& o) const -> Point3d<decltype(x/o.x)> {
      return Point3d<decltype(x/o.x)>(x/o.x, y/o.y, z/o.z);
    }
    auto operator/(const T& o) const -> Point3d<decltype(x/o)> {
      return Point3d<decltype(x/o)>(x/o, y/o, z/o);
    }
  };
}

namespace utils {
  template<typename T>
  fragdock::Point3d<int> round(const fragdock::Point3d<T> &val){
    return fragdock::Point3d<int>(utils::round(val.x), utils::round(val.y), utils::round(val.z));
  }

  /**
  * @param[in] val a point
  * @return a point containing ceiled value
  */
  template<typename T>
  fragdock::Point3d<int> ceili(const fragdock::Point3d<T> &val){
    return fragdock::Point3d<int>(utils::ceili(val.x), utils::ceili(val.y), utils::ceili(val.z));
  }

}

#endif