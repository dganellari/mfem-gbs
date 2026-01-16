#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "mfem.hpp"

namespace mfem {

mfem::real_t u_0(const mfem::Vector &x);
mfem::real_t p_0(const mfem::Vector &x);
void bfield(const mfem::Vector &x, mfem::Vector &v);

} // namespace mfem

#endif // FUNCTIONS_HPP
