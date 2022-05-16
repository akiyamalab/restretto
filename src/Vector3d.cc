#include "Vector3d.hpp"
#include "utils.hpp"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

namespace fragdock {
  void getNormalizeRot(const Vector3d& v1, const Vector3d& _v2, fltype& theta, fltype& phi, fltype& psi) {
    theta = 0;
    phi = 0;
    psi = 0;
    Vector3d vec1 = v1;
    Vector3d vec2 = v1.cross(_v2);
    psi = Vector3d(0, 1, 0).getAngle(Vector3d(vec2.x, vec2.y, 0));
    phi = Vector3d(0, 0, 1).getAngle(Vector3d(0, sqrt(vec2.x * vec2.x + vec2.y * vec2.y), vec2.z));
    if (vec2.x < 0) psi = -psi;
    vec1.rotate(0, phi, psi);
    theta = -Vector3d(1, 0, 0).getAngle(vec1);
    if (vec1.y < 0) theta = -theta;
    // rotate(theta, phi, psi);
  }

  Vector3d getRot(const Vector3d& v1, const Vector3d& _v2) {
    fltype theta, phi, psi;
    getNormalizeRot(v1, _v2, theta, phi, psi);
    return Vector3d(-psi, -phi, -theta);
  }

  void Vector3d::axisRotate(const Vector3d &axis, double th) {
    Vector3d n = axis/axis.abs();
    fltype B = cos(th), C = sin(th);
    fltype A = 1-B;
    fltype ox = x, oy = y, oz = z;
    x = (n.x*n.x*A +     B)*ox + (n.x*n.y*A - n.z*C)*oy + (n.z*n.x*A + n.y*C)*oz;
    y = (n.x*n.y*A + n.z*C)*ox + (n.y*n.y*A +     B)*oy + (n.y*n.z*A - n.x*C)*oz;
    z = (n.z*n.x*A - n.y*C)*ox + (n.y*n.z*A + n.x*C)*oy + (n.z*n.z*A +     B)*oz;
  }

  // zxz rotation: thetaだけz軸, phiだけx軸, psiだけz軸に対して，この順番で回転 (オイラー角で検索)
  void Vector3d::rotate(fltype theta, fltype phi, fltype psi) {
    const fltype cosT = cos(theta), sinT = sin(theta), cosH = cos(phi), sinH = sin(phi), cosS = cos(psi), sinS = sin(psi);
    fltype ox = x, oy = y, oz = z;

    x = (cosT*cosS  - sinT*cosH*sinS) * ox + (-cosT*sinS - sinT*cosH*cosS) * oy + (sinT*sinH ) * oz;
    y = (sinT*cosS  + cosT*cosH*sinS) * ox + (-sinT*sinS + cosT*cosH*cosS) * oy + (-cosT*sinH) * oz;
    z = (sinH*sinS                  ) * ox + (sinH*cosS                  ) * oy + (cosH      ) * oz;
  }
  void Vector3d::rotate(const Vector3d& v) {
    rotate(v.x, v.y, v.z);
  }

  const std::vector<Vector3d> makeRotations60() {
    const fltype phi = (1.0 + sqrt(5.0)) / 2.0;
    const fltype pi = acos(-1.0);

    fragdock::Vector3d poleA(0, phi*phi*phi, phi*phi);
    fragdock::Vector3d poleB(0, 1.0, phi*phi);
    poleA = poleA.unit();
    poleB = poleB.unit();
    fltype p = poleA.getAngle(poleB);
    poleA = Vector3d(1, 0, 0);
    poleB = Vector3d(cos(p), sin(p), 0);

    const int ORDER_NUM = 12;
    const int order[ORDER_NUM] = {0, 0, 1, 1, 1, 2, 1, 2, 2, 2, 3, -1};

    std::vector<Vector3d> rots;
    for(int ord_id = 0; ord_id < ORDER_NUM; ord_id++) {
      int o = order[ord_id];
      for(int i = 0; i < 5; i++) {
        // fltype theta = getAxisAngle(poleA, poleB);
        rots.push_back(getRot(poleA, poleB));
        poleB.axisRotate(poleA, 2.0*pi/5.0);
      }
      poleB.axisRotate(poleA, 2.0*o*pi/5.0);
      poleA.axisRotate(poleB, 4.0*pi/3.0);
    }
    return rots;
  }

  const std::vector<Vector3d> readRotations(const std::string& filename) {
    std::ifstream ifs(filename);
    if (ifs.fail()){
      std::cerr << "opening grid file failed:" << filename << std::endl;
      abort();
    }
    std::vector<Vector3d> rots;
    while (!ifs.eof()) {
      std::string str;
      std::getline(ifs, str);
      std::vector<std::string> vals;
      boost::algorithm::split(vals, str, boost::algorithm::is_any_of(","));
      if ((int)vals.size() < 3) continue;
      boost::algorithm::trim(vals[0]);
      boost::algorithm::trim(vals[1]);
      boost::algorithm::trim(vals[2]);
      rots.push_back(Vector3d(boost::lexical_cast<fltype>(vals[0]),
                              boost::lexical_cast<fltype>(vals[1]),
                              boost::lexical_cast<fltype>(vals[2])));
    }
    return rots;
  }
}
