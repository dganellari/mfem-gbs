# Combined core library
add_library(mfem-gbs core/wave_solver.cpp
                     core/alfven_solver.cpp
                     utils/functions.cpp
)
target_link_libraries(mfem-gbs PUBLIC mfem)

# Make src/ directory available for includes.
# CMAKE_CURRENT_SOURCE_DIR is the absolute path to "./src/", the directory of this CMakeLists file.
# This allows includes like << #include "core/alfven_solver.hpp" >>, from anywhere in the project,
# like in "examples/alfven/alfven_main.cpp".
# Otherwise, we would need to have relative includes, such as:
# include "../../core/alfven_solver.hpp"
target_include_directories(mfem-gbs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Executables
add_executable(wave_app   examples/wave/wave_main.cpp)
add_executable(alfven_app examples/alfven/alfven_main.cpp)

# Link executables
target_link_libraries(wave_app   mfem-gbs)
target_link_libraries(alfven_app mfem-gbs)

# Apply DATA_DIR to all executables
target_compile_definitions(wave_app   PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_compile_definitions(alfven_app PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
