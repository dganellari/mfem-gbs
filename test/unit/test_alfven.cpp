#include <gtest/gtest.h>
#include "mfem.hpp"
#include "core/functions.hpp"
#include <string>

// Test structure:
// TEST(TestSuiteName, TestName) { ...body ...}
// ***  TestSuiteName, TestName should NOT contain underscores 

TEST(AlfvenTest, BasicTest)
{
    EXPECT_EQ(1, 1);
}

// Simple test for initial solution.
TEST(AlfvenTest, InitSolU)
{
    mfem::Vector x(3);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 0.0;

    mfem::real_t val = u_0(x);
    EXPECT_EQ(val, 0.0);
}

// // Simple test for initial solution.
// TEST(AlfvenTest, InitSolP)
// {



// }
