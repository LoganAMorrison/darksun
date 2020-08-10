//
// Created by logan on 8/10/20.
//

#ifndef DARKSUN_STIFF_COMMON_HPP
#define DARKSUN_STIFF_COMMON_HPP

#include <functional>

namespace stiff {

struct StiffLinAlg {
  int mle, mue, mbjac, mbb, mdiag, mdiff, mbdiag;
};

using F_fcn = std::function<void(int *, double *, double *, double *)>;
using F_jac = std::function<void(int *, double *, double *, double *, int *)>;
using F_mas = std::function<void(int *, double *, int *)>;

} // namespace stiff

#endif // DARKSUN_STIFF_COMMON_HPP
