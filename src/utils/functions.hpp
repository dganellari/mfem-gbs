#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include "mfem.hpp"

constexpr double pi = 3.14159265358979323846;

mfem::real_t u_0(const mfem::Vector &x);
mfem::real_t p_0(const mfem::Vector &x);
void bfield(const mfem::Vector &x, mfem::Vector &v);

#endif // FUNCTIONS_HPP
