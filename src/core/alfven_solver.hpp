#ifndef ALFVEN_SOLVER_HPP
#define ALFVEN_SOLVER_HPP

#include "mfem.hpp"

namespace mfem {

class AlfvenOperator
{
protected:
    FiniteElementSpace &fespace_u;
    FiniteElementSpace &fespace_p;

    // Essential boundary dofs (for pressure p)
    Array<int> ess_tdof_p;

    BilinearForm *M;
    BilinearForm *N_par;
    BilinearForm *N_full;
    
    MixedBilinearForm *blf_E;

    SparseMatrix M_mat;
    SparseMatrix N_n_mat;
    SparseMatrix E_mat;
    SparseMatrix *F_mat;  // Transpose of E_mat

    BlockOperator *A;
    Array<int> offsets;

    real_t dt_over_two;

public:
    AlfvenOperator(FiniteElementSpace &fespace_u_,
                   FiniteElementSpace &fespace_p_,
                   real_t dt);

    void AssembleSystem();

    void FormRHS(const Vector &u_old, const Vector &p_old, Vector &b) const;

    BlockOperator& GetSystemOperator() { return *A; }

    Array<int>& GetEssentialTDofs() { return ess_tdof_p; }

    real_t ComputeEnergy(const Vector &u, const Vector &p) const;

    ~AlfvenOperator();
};

} // namespace mfem

#endif