#include "functions.hpp"

mfem::real_t u_0(const mfem::Vector &x) {
    
    double X = x(0);
    double Y = x(1);
    double Z = x(2);
    // return std::sin(pi*X) * std::sin(pi*Y) * std::sin(pi*Z);
    return std::sin(pi*X) * std::sin(pi*Y) * (1.- std::cos(pi*Z));
}

mfem::real_t p_0(const mfem::Vector &x) {
    return 0.;
}

void bfield(const mfem::Vector &x, mfem::Vector &returnvalue) { 
   
    returnvalue(0) = 0.0; //+ 0.0*x(0);
    returnvalue(1) = 0.0; //+ 0.0*x(0);
    returnvalue(2) = 1.0; //+ 0.0*x(0);
}
