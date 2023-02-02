#include "common.hpp"
#include <cmath>
#include <iostream>
#ifndef VECTOR3D_H_
#define VECTOR3D_H_

namespace fragdock {
  namespace {
    fltype ACOS(fltype a) {
      if (a >= 1.0) return acos(1.0);
      if (a <= -1.0) return acos(-1.0);
      return acos(a);
    }
  }

  class Vector3d {
  public:
    fltype x, y, z;
    Vector3d() { x = y = z = 0.0; }
    Vector3d(fltype x, fltype y, fltype z) : x(x), y(y), z(z) {}

    fltype dot(const Vector3d &o) const { return x*o.x + y*o.y + z*o.z; }
    Vector3d cross(const Vector3d &o) const { return Vector3d(y*o.z - z*o.y, z*o.x - x*o.z, x*o.y - y*o.x); }
    fltype abs() const { return sqrt(norm()); }
    fltype norm() const { return dot(*this); }
    void axisRotate(const Vector3d &axis, double th);
    Vector3d unit() { fltype ab = this->abs(); return Vector3d(x/ab, y/ab, z/ab); }
    fltype getAngle(const Vector3d &o) const {
      if (this->norm() < EPS || o.norm() < EPS) {
        // std::cerr << "zero vector in Vector3d::getAngle" << std::endl;
        return 0;
      }
      return ACOS(dot(o)/sqrt(this->norm()*o.norm()));
    }
    void rotate(fltype theta, fltype phi, fltype psi);
    void rotate(const Vector3d& v);
    void translate(const Vector3d &o) { x += o.x, y += o.y, z += o.z; }
    Vector3d operator-(const Vector3d &o) const { return Vector3d(x-o.x, y-o.y, z-o.z); }
    Vector3d operator-=(const Vector3d &o) { x -= o.x; y -= o.y; z -= o.z; return *this; }
    Vector3d operator+(const Vector3d &o) const { return Vector3d(x+o.x, y+o.y, z+o.z); }
    Vector3d operator+=(const Vector3d &o) { x += o.x; y += o.y; z += o.z; return *this; }
    Vector3d operator*(const fltype o) const { return Vector3d(x*o, y*o, z*o); }
    Vector3d operator*=(const fltype o) { x *= o, y *= o, z *= o; return *this; }
    Vector3d operator/(const fltype o) const { return *this*(1.0/o); }
    Vector3d operator/=(const fltype o) { x /= o, y /= o, z /= o; return *this; }
    Vector3d operator-() const { return Vector3d(-x, -y, -z); }
    void print() const { std::cout << "(" << x << ", " << y << ", " << z << ")" << std::endl; }
    friend std::ostream& operator<< (std::ostream& os, const Vector3d& v) { os << "(" << v.x << ", " << v.y << ", " << v.z << ")" << std::endl; return os; }
  };

  const std::vector<Vector3d> makeRotations60();
  const std::vector<Vector3d> readRotations(const std::string& rotangs_file);
}


#endif