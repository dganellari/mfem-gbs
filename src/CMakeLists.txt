# Combined core library
add_library(mfem-gbs core/wave_solver.cpp
                     core/alfven_solver.cpp
                     utils/functions.cpp
)
target_link_libraries(mfem-gbs PUBLIC mfem)

# Executables
add_executable(wave_app   examples/wave/wave_main.cpp)
add_executable(alfven_app examples/alfven/alfven_main.cpp)

# Link executables
target_link_libraries(wave_app   mfem-gbs)
target_link_libraries(alfven_app mfem-gbs)

# Apply DATA_DIR to all executables
target_compile_definitions(wave_app   PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
target_compile_definitions(alfven_app PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
