// MFEM Wave Mini-App
// Simplified main using modular wave solver

#include "mfem.hpp"
#include "wave_solver.hpp"
#include <iostream>

using namespace std;
using namespace mfem;

int main()
{
   // Parameters
   const char *mesh_file = "star.mesh";
   int ref_levels = 1;
   int order = 2;
   int ode_solver_type = 10;
   real_t t_final = 0.1;
   real_t dt = 1.0e-2;
   real_t speed = 1.0;
   bool dirichlet = true;

   cout << "MFEM Wave Mini-App" << endl;

   // Read mesh
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // ODE solver
   SecondOrderODESolver *ode_solver = SecondOrderODESolver::Select(ode_solver_type);

   // Refine mesh
   for (int lev = 0; lev < ref_levels; lev++)
   {
      mesh.UniformRefinement();
   }

   // FE space
   H1_FECollection fe_coll(order, dim);
   FiniteElementSpace fespace(&mesh, &fe_coll);

   cout << "Number of unknowns: " << fespace.GetTrueVSize() << endl;

   GridFunction u_gf(&fespace);
   GridFunction dudt_gf(&fespace);

   // Initial conditions
   FunctionCoefficient u_0(InitialSolution);
   u_gf.ProjectCoefficient(u_0);
   Vector u;
   u_gf.GetTrueDofs(u);

   FunctionCoefficient dudt_0(InitialRate);
   dudt_gf.ProjectCoefficient(dudt_0);
   Vector dudt;
   dudt_gf.GetTrueDofs(dudt);

   // Wave operator
   Array<int> ess_bdr(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0);
   if (dirichlet) ess_bdr = 1; else ess_bdr = 0;
   WaveOperator oper(fespace, ess_bdr, speed);

   // Time integration
   ode_solver->Init(oper);
   real_t t = 0.0;

   int steps = t_final / dt;
   for (int ti = 1; ti <= steps; ti++)
   {
      ode_solver->Step(u, dudt, t, dt);
      cout << "Step " << ti << ", t = " << t << endl;
      oper.SetParameters(u);
   }

   // Save final solution
   u_gf.SetFromTrueDofs(u);
   ofstream osol("final_solution.gf");
   u_gf.Save(osol);

   cout << "Simulation complete. Solution saved to final_solution.gf" << endl;

   delete ode_solver;

   return 0;
}