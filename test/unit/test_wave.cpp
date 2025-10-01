#include <gtest/gtest.h>
#include "mfem.hpp"
#include "core/wave_solver.hpp"
#include <string>

using namespace mfem;

// Test fixture for common setup
class WaveTestFixture : public ::testing::Test {
protected:
    std::string mesh_file;
    Mesh* mesh;
    H1_FECollection* fe_coll;
    FiniteElementSpace* fespace;
    
    void SetUp() override {
        mesh_file = std::string(DATA_DIR) + "/star.mesh";
        mesh = new Mesh(mesh_file.c_str(), 1, 1);
        fe_coll = new H1_FECollection(1, mesh->Dimension());
        fespace = new FiniteElementSpace(mesh, fe_coll);
    }
    
    void TearDown() override {
        delete fespace;
        delete fe_coll;
        delete mesh;
    }
};

// Test initial solution
TEST(WaveTest, InitialSolution)
{
   Vector x(2);
   x[0] = 0.0; x[1] = 0.0;
   real_t val = InitialSolution(x);
   EXPECT_GT(val, 0.0);
   
   // Test at different points
   x[0] = 1.0; x[1] = 1.0;
   real_t val2 = InitialSolution(x);
   EXPECT_GE(val2, 0.0);
}

// Test initial rate
TEST(WaveTest, InitialRate)
{
   Vector x(2);
   x[0] = 0.0; x[1] = 0.0;
   real_t val = InitialRate(x);
   EXPECT_EQ(val, 0.0);
}

// Test operator creation using fixture
TEST_F(WaveTestFixture, OperatorCreation)
{
   Array<int> ess_bdr(mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0);
   ess_bdr = 1;

   WaveOperator oper(*fespace, ess_bdr, 1.0);
   EXPECT_EQ(oper.Height(), fespace->GetTrueVSize());
   EXPECT_EQ(oper.Width(), fespace->GetTrueVSize());
}

// Test a single time step
TEST_F(WaveTestFixture, TimeStep)
{
   GridFunction u_gf(fespace);
   FunctionCoefficient u_0(InitialSolution);
   u_gf.ProjectCoefficient(u_0);
   Vector u;
   u_gf.GetTrueDofs(u);

   Vector dudt(u.Size());
   dudt = 0.0;

   Array<int> ess_bdr(mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0);
   ess_bdr = 1;
   WaveOperator oper(*fespace, ess_bdr, 1.0);

   SecondOrderODESolver *ode_solver = SecondOrderODESolver::Select(10);
   ASSERT_NE(ode_solver, nullptr) << "Failed to create ODE solver";
   
   ode_solver->Init(oper);

   real_t t = 0.0;
   real_t dt = 1e-2;
   Vector u_initial = u;  // Store initial state
   
   ode_solver->Step(u, dudt, t, dt);

   // Solution should change after time step
   EXPECT_NE(u[0], u_initial[0]);

   delete ode_solver;
}

// Test with different speed
TEST_F(WaveTestFixture, DifferentSpeed)
{
   Array<int> ess_bdr(mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0);
   ess_bdr = 1;

   WaveOperator oper1(*fespace, ess_bdr, 1.0);
   WaveOperator oper2(*fespace, ess_bdr, 2.0);
   
   // Both operators should have same size
   EXPECT_EQ(oper1.Height(), oper2.Height());
   EXPECT_EQ(oper1.Width(), oper2.Width());
}

// Parameterized test example for future expansion
class WaveSpeedTest : public WaveTestFixture,
                      public ::testing::WithParamInterface<double> {
};

TEST_P(WaveSpeedTest, VariousSpeed)
{
   double speed = GetParam();
   Array<int> ess_bdr(mesh->bdr_attributes.Size() ? mesh->bdr_attributes.Max() : 0);
   ess_bdr = 1;
   
   WaveOperator oper(*fespace, ess_bdr, speed);
   EXPECT_EQ(oper.Height(), fespace->GetTrueVSize());
}

INSTANTIATE_TEST_SUITE_P(SpeedRange, WaveSpeedTest,
                        ::testing::Values(0.5, 1.0, 2.0, 5.0));