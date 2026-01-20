#include <gtest/gtest.h>
#include "mfem.hpp"
#include "core/alfven_solver.hpp"
#include "utils/functions.hpp"
#include <string>

// Test structure:
// TEST(TestSuiteName, TestName) { ...body ...}
// ***  TestSuiteName, TestName should NOT contain underscores 

// TEST   is for simple tests without shared setup/teardown
// TEST_F is for tests that share common setup/teardown via a fixture class

// Simple test for initial solution.
TEST(AlfvenTest, InitialSolutionU)
{
    mfem::Vector x(3);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;
    mfem::real_t val1 = mfem::u_0(x);
    EXPECT_EQ(val1, 0.0);

    // Test at another (non-zero) point
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 0.5;
    mfem::real_t val2 = mfem::u_0(x);
    EXPECT_GT(val2, 0.0); // should be non-zero
}

// Simple test for initial solution.
TEST(AlfvenTest, InitialSolutionP)
{
    mfem::Vector x(3);
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 0.5;
    mfem::real_t val = mfem::p_0(x); // p_0 is zero everywhere
    EXPECT_EQ(val, 0.0) << "Initial pressure should be zero everywhere";
}

// Test bfield
TEST(AlfvenTest, MagneticField)
{
    mfem::Vector x(3);
    x[0] = 0.5;
    x[1] = 0.5;
    x[2] = 0.5;
    mfem::Vector b(3);
    mfem::bfield(x, b);
    
    // to revise expected values based on bfield definition
    EXPECT_EQ(b[0], 0.0);
    EXPECT_EQ(b[1], 0.0);
    EXPECT_EQ(b[2], 1.0);
}

// Common setup (mesh, FE space) for fixture tests
class AlfvenTestFixture : public ::testing::Test {
protected:
    // Static mesh shared across all tests (loaded once)
    static mfem::Mesh* mesh;
    
    // Per-test FE collections and spaces
    mfem::H1_FECollection* fe_coll_u;
    mfem::H1_FECollection* fe_coll_p; // could use the same FE for both u,p
    mfem::FiniteElementSpace* fespace_u;
    mfem::FiniteElementSpace* fespace_p;
    
    // Called once before all tests in this suite
    static void SetUpTestSuite() {
        std::string mesh_file = std::string(DATA_DIR) + "/ref-cube.mesh";
        // Load mesh once for all tests
        mesh = new mfem::Mesh(mesh_file.c_str(), 1, 1);
    }
    
    // Called once after all tests in this suite
    static void TearDownTestSuite() {
        delete mesh;
    }
    
    // SetUp() is called before EACH test.
    // In principle, I could move everthing to SetUpTestSuite(),
    //  but keeping FE spaces separate per test is cleaner
    void SetUp() override {
        int dim = mesh->Dimension();
        
        // H1 space for velocity u
        fe_coll_u = new mfem::H1_FECollection(1, dim);
        fespace_u = new mfem::FiniteElementSpace(mesh, fe_coll_u);
        
        // H1 space for pressure p
        fe_coll_p = new mfem::H1_FECollection(1, dim);
        fespace_p = new mfem::FiniteElementSpace(mesh, fe_coll_p);
    }
    
    // Called after each test
    void TearDown() override {
        delete fe_coll_u;
        delete fe_coll_p;
        delete fespace_u;
        delete fespace_p;
    }
};

// Initialize static mesh
mfem::Mesh* AlfvenTestFixture::mesh = nullptr;

// Test: AlfvenOperator constructs successfully with correct block dimensions
TEST_F(AlfvenTestFixture, OperatorCreation)
{
    mfem::real_t dt = 0.01;
    mfem::AlfvenOperator oper(*fespace_u, *fespace_p, dt);
    
    mfem::BlockOperator& A = oper.GetSystemOperator();
    
    // Check dimensions
    int expected_size = fespace_u->GetNDofs() + fespace_p->GetNDofs();
    EXPECT_EQ(A.Height(), expected_size);
    EXPECT_EQ(A.Width(), expected_size);
}

// Test: Essential boundary DOFs are correctly identified and within valid range
TEST_F(AlfvenTestFixture, EssentialDOFs)
{
    mfem::real_t dt = 0.01;
    mfem::AlfvenOperator oper(*fespace_u, *fespace_p, dt);
    
    mfem::Array<int>& ess_tdofs = oper.GetEssentialTDofs();
    
    // Essential DOFs should exist
    // boundary conditions for pressure, p is in the 2nd block
    EXPECT_GE(ess_tdofs.Size(), 0);
    
    // All DOF indices should be valid
    int total_size = fespace_u->GetNDofs() + fespace_p->GetNDofs();
    for (int i = 0; i < ess_tdofs.Size(); i++)
    {
        EXPECT_GT(ess_tdofs[i], 0);
        EXPECT_LT(ess_tdofs[i], total_size);
        // std::cout << "Essential DOF " << i << ": " << ess_tdofs[i] << std::endl;
    }
}

// Test: FormRHS produces non-zero output for non-zero input
TEST_F(AlfvenTestFixture, RHSFormation)
{
    mfem::real_t dt = 0.01;
    mfem::AlfvenOperator oper(*fespace_u, *fespace_p, dt);

    int u_size = fespace_u->GetNDofs();
    int p_size = fespace_p->GetNDofs();
    int total_size = u_size + p_size;

    mfem::Vector u(u_size);
    mfem::Vector p(p_size);

    u = 1.0;
    p = 0.5;
    
    // Compute RHS for baseline state
    mfem::Vector b1(total_size);
    oper.FormRHS(u, p, b1);
    mfem::real_t norm1 = b1.Norml2();
    EXPECT_GT(norm1, 0.0) << "RHS should be non-zero for non-zero input";

    // Double the input state
    u *= 2.0;
    p *= 2.0;

    mfem::Vector b2(total_size);
    oper.FormRHS(u, p, b2);
    mfem::real_t norm2 = b2.Norml2();

    // Doubling inputs should double outputs (linearity property)
    mfem::real_t ratio = norm2 / norm1;
    EXPECT_NEAR(ratio, 2.0, 1e-10) << "FormRHS must be linear operator";
}

// Test: Energy should remain bounded and scale predictably with state magnitude
TEST_F(AlfvenTestFixture, EnergyComputation)
{
    mfem::real_t dt = 0.01;
    mfem::AlfvenOperator oper(*fespace_u, *fespace_p, dt);

    int u_size = fespace_u->GetNDofs();
    int p_size = fespace_p->GetNDofs();

    // Test energy at different state magnitudes (including zero)
    std::vector<mfem::real_t> energies;
    
    // Allocate vectors once
    mfem::Vector u(u_size);
    mfem::Vector p(p_size);
    
    for (int scale = 0; scale <= 10; ++scale) {
        u = 0.10 * scale;
        p = 0.05 * scale;
        
        mfem::real_t energy = oper.ComputeEnergy(u, p);
        
        // Zero state should have zero energy
        if (scale == 0) {
            // Verifies that the two DOUBLE values are approximately equal
            EXPECT_DOUBLE_EQ(energy, 0.0);
        }
        
        // Energy must be finite (no NaN/Inf)
        EXPECT_FALSE(std::isnan(energy)) << "Energy should not be NaN at scale " << scale;
        EXPECT_FALSE(std::isinf(energy)) << "Energy should not be infinite at scale " << scale;
        
        // Energy should increase with state magnitude (quadratic form property)
        if (scale > 0) {
            EXPECT_GT(energy, energies[scale-1]) 
                << "Energy should grow with state amplitude";
        }
        energies.push_back(energy);
    }
}

// Test: Verify that E and F matrices are transposes of each other
TEST_F(AlfvenTestFixture, CouplingBlocks)
{
    mfem::real_t dt = 0.01;
    mfem::AlfvenOperator oper(*fespace_u, *fespace_p, dt);
    
    const mfem::SparseMatrix &E = oper.Get_E_Matrix();
    const mfem::SparseMatrix &F = oper.Get_F_Matrix();

    // Verify F = E^T by checking dimensions
    EXPECT_EQ(E.Height(), F.Width());
    EXPECT_EQ(E.Width(),  F.Height());

    // Verify F = E^T numerically: F*x should equal E^T*x
    mfem::Vector x(E.Width());
    x.Randomize(1);

    mfem::Vector y_from_F(F.Height());           // y = F*x
    mfem::Vector y_from_E_transpose(E.Height()); // y = E^T*x

    F.Mult(x, y_from_F);                    // y = F*x
    E.MultTranspose(x, y_from_E_transpose); // y = E^T*x

    y_from_F -= y_from_E_transpose;  // Should be zero
    EXPECT_LT(y_from_F.Norml2(), 1e-12) << "F should equal E^T";
}

// Test: Different time steps should affect RHS differently
// The coupling (off-diagonal) matrices E, F (F = E^T) depend on dt (scaled by dt/2)
TEST_F(AlfvenTestFixture, TimestepScaling)
{
    int u_size = fespace_u->GetNDofs();
    int p_size = fespace_p->GetNDofs();
    int total_size = u_size + p_size;

    mfem::Vector u(u_size);
    mfem::Vector p(p_size);

    u = 1.0;
    p = 1.0;

    const mfem::real_t small_dt = 0.001;
    const mfem::real_t large_dt = 0.100;    // 100x larger

    // Small time step
    mfem::AlfvenOperator oper1(*fespace_u, *fespace_p, small_dt);
    mfem::Vector b1(total_size);
    oper1.FormRHS(u, p, b1);

    // Large time step
    mfem::AlfvenOperator oper2(*fespace_u, *fespace_p, large_dt);
    mfem::Vector b2(total_size);
    oper2.FormRHS(u, p, b2);

    // RHS should differ because E_mat *= -dt/2 in AssembleSystem()
    mfem::real_t norm1 = b1.Norml2();
    mfem::real_t norm2 = b2.Norml2();

    EXPECT_NE(norm1, norm2) << "Different dt should produce different RHS";
    // Larger dt should give larger coupling contribution
    EXPECT_GT(norm2, norm1) << "Larger dt should increase coupling strength";
}
