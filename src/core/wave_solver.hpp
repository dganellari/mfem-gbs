#ifndef WAVE_SOLVER_HPP
#define WAVE_SOLVER_HPP

#include "mfem.hpp"

namespace mfem {

class WaveOperator : public SecondOrderTimeDependentOperator
{
protected:
   FiniteElementSpace &fespace;
   Array<int> ess_tdof_list;

   BilinearForm *M;
   BilinearForm *K;

   SparseMatrix Mmat, Kmat;
   SparseMatrix *T;
   real_t current_dt;

   CGSolver M_solver;
   DSmoother M_prec;

   CGSolver T_solver;
   DSmoother T_prec;

   Coefficient *c2;
   mutable Vector z;

public:
   WaveOperator(FiniteElementSpace &f, Array<int> &ess_bdr, real_t speed);

   using SecondOrderTimeDependentOperator::Mult;
   void Mult(const Vector &u, const Vector &du_dt,
             Vector &d2udt2) const override;

   using SecondOrderTimeDependentOperator::ImplicitSolve;
   void ImplicitSolve(const real_t fac0, const real_t fac1,
                      const Vector &u, const Vector &dudt, Vector &d2udt2) override;

   void SetParameters(const Vector &u);

   ~WaveOperator() override;
};

real_t InitialSolution(const Vector &x);
real_t InitialRate(const Vector &x);

} // namespace mfem

#endif