#include "functions.hpp"

namespace mfem {

mfem::real_t u_0(const mfem::Vector &x) {
    
    double X = x(0);
    double Y = x(1);
    double Z = x(2);
    return std::sin(M_PI*X) * std::sin(M_PI*Y) * (1.- std::cos(M_PI*Z));
    // M_PI is provided by <cmath>
}

mfem::real_t p_0(const mfem::Vector &x) {
    return 0.;
}

void bfield(const mfem::Vector &x, mfem::Vector &returnvalue) { 
   
    returnvalue(0) = 0.0; //+ 0.0*x(0);
    returnvalue(1) = 0.0; //+ 0.0*x(0);
    returnvalue(2) = 1.0; //+ 0.0*x(0);
}

} // namespace mfem
