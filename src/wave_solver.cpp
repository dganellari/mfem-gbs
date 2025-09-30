#include "wave_solver.hpp"

namespace mfem {

WaveOperator::WaveOperator(FiniteElementSpace &f,
                           Array<int> &ess_bdr, real_t speed)
   : SecondOrderTimeDependentOperator(f.GetTrueVSize(), (real_t) 0.0),
     fespace(f), M(NULL), K(NULL), T(NULL), current_dt(0.0), z(height)
{
   c2 = new ConstantCoefficient(speed*speed);
   K = new BilinearForm(&fespace);
   K->AddDomainIntegrator(new DiffusionIntegrator(*c2));
   K->Assemble();

   M = new BilinearForm(&fespace);
   M->AddDomainIntegrator(new MassIntegrator());
   M->Assemble();

   fespace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);
   K->FormSystemMatrix(ess_tdof_list, Kmat);
   M->FormSystemMatrix(ess_tdof_list, Mmat);

   const real_t rel_tol = 1e-8;
   M_solver.iterative_mode = false;
   M_solver.SetRelTol(rel_tol);
   M_solver.SetAbsTol(0.0);
   M_solver.SetMaxIter(30);
   M_solver.SetPrintLevel(0);
   M_solver.SetPreconditioner(M_prec);
   M_solver.SetOperator(Mmat);

   T_solver.iterative_mode = false;
   T_solver.SetRelTol(rel_tol);
   T_solver.SetAbsTol(0.0);
   T_solver.SetMaxIter(100);
   T_solver.SetPrintLevel(0);
   T_solver.SetPreconditioner(T_prec);
}

void WaveOperator::Mult(const Vector &u, const Vector &du_dt,
                        Vector &d2udt2) const
{
   K->FullMult(u, z);
   z.Neg();
   z.SetSubVector(ess_tdof_list, 0.0);
   M_solver.Mult(z, d2udt2);
   d2udt2.SetSubVector(ess_tdof_list, 0.0);
}

void WaveOperator::ImplicitSolve(const real_t fac0, const real_t fac1,
                                 const Vector &u, const Vector &dudt, Vector &d2udt2)
{
   if (!T)
   {
      T = Add(1.0, Mmat, fac0, Kmat);
      T_solver.SetOperator(*T);
   }
   K->FullMult(u, z);
   z.Neg();
   z.SetSubVector(ess_tdof_list, 0.0);
   T_solver.Mult(z, d2udt2);
   d2udt2.SetSubVector(ess_tdof_list, 0.0);
}

void WaveOperator::SetParameters(const Vector &u)
{
   delete T;
   T = NULL;
}

WaveOperator::~WaveOperator()
{
   delete T;
   delete M;
   delete K;
   delete c2;
}

real_t InitialSolution(const Vector &x)
{
   return exp(-x.Norml2()*x.Norml2()*30);
}

real_t InitialRate(const Vector &x)
{
   return 0.0;
}

} // namespace mfem