# MFEM Wave Mini-App

A module MFEM-based solver demonstrating wave equation simulation with CMake and GTests.
Modular design for easy expansion, with Test-Drive Development (TDD) workflow.

## Features
- **Modular Architecture**: Separates core solver logic from application
- **TDD Workflow**: Comprehensive test suite with GTest
- **Flexible Configuration**: Command-line argument parsing
- **Scalable Design**: Ready for expansion to complex simulations

## Quick Start

### Prerequisites
```bash
# Install MFEM (via Spack recommended)
spack install mfem

# Or build MFEM from source
git clone https://github.com/mfem/mfem.git
cd mfem && make serial -j
```

### Build
```bash
mkdir build && cd build
cmake .. -DCMAKE_PREFIX_PATH=/path/to/mfem/install
make -j
```

### Run
```bash
# Run with defaults
./wave_app

# Custom parameters
./wave_app -m ../data/star.mesh -r 2 -o 3 -tf 1.0 -dt 0.001

# Run tests
ctest --output-on-failure
```

## TDD Workflow
1. 🔴 **Red**: Write a failing test in `test/unit/test_wave.cpp`
2. 🟢 **Green**: Implement minimal code to pass
3. 🔵 **Refactor**: Improve code while keeping tests green
4. 🔄 **Repeat**: Run `ctest` continuously during development


## Project Structure
```
├── src/              # Source code
│   ├── core/         # Core solver logic
│   │   ├── wave_solver.hpp
│   │   └── wave_solver.cpp
│   ├── solvers/      # Additional solvers (future)
│   ├── utils/        # Utility functions (future)
│   └── main.cpp      # Application entry point
├── test/             # Tests
│   ├── unit/         # Unit tests
│   │   └── test_wave.cpp
│   └── integration/  # Integration tests (future)
├── data/             # Mesh files and data
│   └── star.mesh
└── docs/             # Documentation (future)
```

## Future Enhancements
- [ ] Parallel (MPI) support
- [ ] GPU acceleration (CUDA/HIP)
- [ ] Visualization integration
- [ ] Performance benchmarks
- [ ] CI/CD pipeline

## Contributing
Follow TDD practices: tests first, then implementation.

## License
[Your License]