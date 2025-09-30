add_executable(test_wave test_wave.cpp)
target_link_libraries(test_wave GTest::gtest_main wave_solver)
target_compile_definitions(test_wave PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")

add_test(NAME wave_test COMMAND test_wave)
include(GoogleTest)
gtest_discover_tests(test_wave TEST_LIST test_list)
set_tests_properties(${test_list} PROPERTIES ENVIRONMENT "MFEM_DATA_DIR=${CMAKE_SOURCE_DIR}/data")