#ifndef DARKSUN_STIFF_COMMON_HPP
#define DARKSUN_STIFF_COMMON_HPP

#include "darksun/stiff/vector_matrix.hpp"
#include <functional>

namespace darksun::stiff {
template <int N>
using OdeFun =
    std::function<void(double, const Vector<double, N> &, Vector<double, N> &)>;
template <int N>
using OdeJac = std::function<void(double, const Vector<double, N> &,
                                  Matrix<double, N, N> &)>;

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_COMMON_HPP
