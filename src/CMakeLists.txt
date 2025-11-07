add_library(wave_solver core/wave_solver.cpp)
target_link_libraries(wave_solver PUBLIC mfem)

add_executable(wave_app main.cpp)
target_link_libraries(wave_app wave_solver)
target_compile_definitions(wave_app PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")

#################################################################################
add_library(functions core/functions.cpp)
target_link_libraries(functions PUBLIC mfem)

add_executable(alfven alfven.cpp)
target_link_libraries(alfven functions)
target_compile_definitions(alfven PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")
