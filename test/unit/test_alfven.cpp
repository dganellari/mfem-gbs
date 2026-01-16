#include <gtest/gtest.h>
#include "mfem.hpp"
#include "core/alfven_solver.hpp"
#include "utils/functions.hpp"
#include <string>

// Test structure:
// TEST(TestSuiteName, TestName) { ...body ...}
// ***  TestSuiteName, TestName should NOT contain underscores 

TEST(AlfvenTest, BasicTest)
{
    EXPECT_EQ(1, 1);
}

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
