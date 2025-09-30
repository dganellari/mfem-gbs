#include <gtest/gtest.h>
#include "mfem.hpp"
#include "../src/wave_solver.hpp"

using namespace mfem;

// Test operator creation
TEST(WaveTest, OperatorCreation)
{
   Mesh mesh("star.mesh", 1, 1);
   H1_FECollection fe_coll(1, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fe_coll);
   Array<int> ess_bdr(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0);
   ess_bdr = 1;

   WaveOperator oper(fespace, ess_bdr, 1.0);
   EXPECT_EQ(oper.Height(), fespace.GetTrueVSize());
}

// Test initial solution
TEST(WaveTest, InitialSolution)
{
   Vector x(2);
   x[0] = 0.0; x[1] = 0.0;
   real_t val = InitialSolution(x);
   EXPECT_GT(val, 0.0);
}

// Test initial rate
TEST(WaveTest, InitialRate)
{
   Vector x(2);
   x[0] = 0.0; x[1] = 0.0;
   real_t val = InitialRate(x);
   EXPECT_EQ(val, 0.0);
}

// Test a single time step
TEST(WaveTest, TimeStep)
{
   Mesh mesh("star.mesh", 1, 1);
   H1_FECollection fe_coll(1, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fe_coll);

   GridFunction u_gf(&fespace);
   FunctionCoefficient u_0(InitialSolution);
   u_gf.ProjectCoefficient(u_0);
   Vector u;
   u_gf.GetTrueDofs(u);

   Vector dudt(u.Size());
   dudt = 0.0;

   Array<int> ess_bdr(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0);
   ess_bdr = 1;
   WaveOperator oper(fespace, ess_bdr, 1.0);

   SecondOrderODESolver *ode_solver = SecondOrderODESolver::Select(10);
   ode_solver->Init(oper);

   real_t t = 0.0;
   real_t dt = 1e-2;
   ode_solver->Step(u, dudt, t, dt);

   EXPECT_NE(u[0], 0.0);

   delete ode_solver;
}

// Test with different speed
TEST(WaveTest, DifferentSpeed)
{
   Mesh mesh("star.mesh", 1, 1);
   H1_FECollection fe_coll(1, mesh.Dimension());
   FiniteElementSpace fespace(&mesh, &fe_coll);
   Array<int> ess_bdr(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0);
   ess_bdr = 1;

   WaveOperator oper1(fespace, ess_bdr, 1.0);
   WaveOperator oper2(fespace, ess_bdr, 2.0);

   // Different speeds should create different operators (height same, but internals differ)
   EXPECT_EQ(oper1.Height(), oper2.Height());
}