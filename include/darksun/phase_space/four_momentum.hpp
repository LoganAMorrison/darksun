#ifndef DARK_SUN_PHASE_SPACE_FOUR_MOMENTUM_HPP
#define DARK_SUN_PHASE_SPACE_FOUR_MOMENTUM_HPP

#include <cmath>
#include <iostream>

namespace darksun {

class FourMomentum {
public:
  double e;
  double px;
  double py;
  double pz;

  double mass();
  double dot(const FourMomentum &);
  FourMomentum operator+(const FourMomentum &);
  FourMomentum operator-(const FourMomentum &);
};

FourMomentum FourMomentum::operator+(const FourMomentum &fv) {
  return FourMomentum{e + fv.e, px + fv.px, py + fv.py, pz + fv.pz};
}

FourMomentum FourMomentum::operator-(const FourMomentum &fv) {
  return FourMomentum{e - fv.e, px - fv.px, py - fv.py, pz - fv.pz};
}

FourMomentum operator+(const FourMomentum &fv1, const FourMomentum &fv2) {
  return FourMomentum{fv1.e + fv2.e, fv1.px + fv2.px, fv1.py + fv2.py,
                      fv1.pz + fv2.pz};
}

FourMomentum operator-(const FourMomentum &fv1, const FourMomentum &fv2) {
  return FourMomentum{fv1.e - fv2.e, fv1.px - fv2.px, fv1.py - fv2.py,
                      fv1.pz - fv2.pz};
}

std::ostream &operator<<(std::ostream &os, const FourMomentum &fv) {
  os << "FourMomentum(" << fv.e << ", " << fv.px << ", " << fv.py << ", "
     << fv.pz << ")";
  return os;
}

double FourMomentum::mass() {
  double m = e * e - px * px - py * py - pz * pz;
  return m >= 0 ? std::sqrt(m) : -std::sqrt(-m);
}

double FourMomentum::dot(const FourMomentum &fv) {
  return e * fv.e - px * fv.px - py * fv.py - pz * fv.pz;
}

double scalar_product(const FourMomentum &fv1, const FourMomentum &fv2) {
  return fv1.e * fv2.e - fv1.px * fv2.px - fv1.py * fv2.py - fv1.pz * fv2.pz;
}

} // namespace darksun

#endif // DARK_SUN_PHASE_SPACE_FOUR_MOMENTUM_HPP
