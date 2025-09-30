cmake_minimum_required(VERSION 3.12)

include(FetchContent)

project(mfem-wave-miniapp CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find MFEM
find_package(MFEM REQUIRED)
message(STATUS "Found MFEM version: ${MFEM_VERSION}")

# Find GTest
find_package(GTest QUIET)
if(GTest_FOUND)
    message(STATUS "GTest found: ${GTest_VERSION}")
else()
    message(STATUS "GTest not found, fetching from source")
endif()

if(NOT GTest_FOUND)
    FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(googletest)
endif()

include(GoogleTest)

enable_testing()

include_directories($PROJECT_SOURCE_DIR}/src)

# Add subdirectories for modularity
add_subdirectory(src)
add_subdirectory(test)

# Enable testing

# If MPI is available, add MPI test (optional)
# if (MFEM_USE_MPI)
#   add_test(NAME wave_test_mpi_np=2
#     COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2
#     $<TARGET_FILE:test_wave>)
# endif()