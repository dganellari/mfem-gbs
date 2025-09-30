cmake_minimum_required(VERSION 3.12)
project(mfem-wave-miniapp CXX)

# Find MFEM
find_package(MFEM REQUIRED)

# Find GTest
find_package(GTest REQUIRED)

# Add executable
add_executable(wave_app src/main.cpp src/wave_solver.cpp)
target_link_libraries(wave_app MFEM::mfem)

# Enable testing
enable_testing()

# Add test executable
add_executable(test_wave test/test_wave.cpp src/wave_solver.cpp)
target_link_libraries(test_wave GTest::gtest_main MFEM::mfem)

# Add test with options (similar to MFEM examples)
add_test(NAME wave_test COMMAND test_wave)

# If MPI is available, add MPI test (optional)
if (MFEM_USE_MPI)
  add_test(NAME wave_test_mpi_np=2
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
    $<TARGET_FILE:test_wave>)
endif()