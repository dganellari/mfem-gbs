# MFEM Wave Mini-App

A standalone mini-application based on MFEM Example 23, demonstrating wave equation simulation with CMake and GTests. Modular design for easy expansion, with TDD support.

## Structure
- `src/wave_solver.hpp/cpp`: Core wave simulation logic (WaveOperator class).
- `src/main.cpp`: Simple main function using the wave solver.
- `test/test_wave.cpp`: Unit/integration tests for the wave solver.

## TDD Workflow
This project follows Test Driven Development:
1. Write failing tests in `test/test_wave.cpp` for new features.
2. Implement minimal code in `src/` to pass tests.
3. Refactor while keeping tests passing.
4. Run `./run_tests.sh` frequently for feedback.

## Requirements
- MFEM library installed (via Spack, CMake, etc.)
- GTest installed
- CMake 3.12+

## Build Instructions
1. Ensure MFEM and GTest are installed and findable by CMake.
2. Run:
   ```
   mkdir build && cd build
   cmake ..
   make
   ```
3. Run the app: `./wave_app`
4. Run tests: `ctest` or `./test_wave` or `./run_tests.sh`

## Learned from MFEM Examples CMakeLists
- **Conditional Tests**: Added MPI test if MFEM_USE_MPI is enabled (like MFEM's parallel examples).
- **Test Options**: Tests can include command-line flags for different scenarios (e.g., device backends).
- **Modular Builds**: Subdirectories for features (e.g., if adding PETSc-specific tests later).
- **Skip Logic**: Skip tests based on features (e.g., if single precision, skip certain tests).

For expansion, consider adding device tests (CUDA/HIP) or solver-specific variants.

## Notes
- Assumes `star.mesh` is in the run directory (copy from MFEM data/).
- Designed for expansion: add more solvers, options, etc.