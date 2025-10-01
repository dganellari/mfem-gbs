cmake_minimum_required(VERSION 3.12)

include(FetchContent)

project(mfem-wave-miniapp 
    VERSION 0.1.0
    DESCRIPTION "MFEM Wave Equation Miniapp with TDD"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build options
option(BUILD_TESTING "Build tests" ON)
# option(BUILD_EXAMPLES "Build examples" ON)

# Find MFEM
find_package(MFEM REQUIRED)
message(STATUS "Found MFEM version: ${MFEM_VERSION}")
message(STATUS "MFEM_USE_MPI: ${MFEM_USE_MPI}")

# Find GTest (only if testing is enabled)
if(BUILD_TESTING)
    find_package(GTest QUIET)
    if(GTest_FOUND)
        message(STATUS "GTest found: ${GTest_VERSION}")
    else()
        message(STATUS "GTest not found, fetching from source")
        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
            DOWNLOAD_EXTRACT_TIMESTAMP TRUE
        )
        # Don't install GTest
        set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
        set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable(googletest)
    endif()
    
    include(GoogleTest)
    enable_testing()
endif()

# Include directories
include_directories(${PROJECT_SOURCE_DIR}/src)

# Add subdirectories
add_subdirectory(src)

if(BUILD_TESTING)
    add_subdirectory(test)
endif()

# Installation rules
install(TARGETS wave_app DESTINATION bin)