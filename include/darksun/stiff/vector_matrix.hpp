#ifndef DARKSUN_STIFF_VECTOR_MATRIX_HPP
#define DARKSUN_STIFF_VECTOR_MATRIX_HPP

#include <Eigen/Dense>

namespace darksun ::stiff {

template <typename T, int SizeN = Eigen::Dynamic, int SizeM = Eigen::Dynamic>
using Matrix = Eigen::Matrix<T, SizeN, SizeM>;

template <typename T, int SizeN = Eigen::Dynamic>
using Vector = Eigen::Matrix<T, SizeN, 1>;

} // namespace darksun::stiff

#endif // DARKSUN_STIFF_VECTOR_MATRIX_HPP
