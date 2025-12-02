add_library(wave_solver core/wave_solver.cpp)
target_link_libraries(wave_solver PUBLIC mfem)

# Original Alfven app
add_executable(wave_app alfven_dec.cpp)
target_link_libraries(wave_app wave_solver)
target_compile_definitions(wave_app PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# CUDA-enabled Alfven app
add_executable(alfven_cuda_app alfven_dec_cuda.cpp)
target_link_libraries(alfven_cuda_app wave_solver)
target_compile_definitions(alfven_cuda_app PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
