# wave app
add_executable(test_wave unit/test_wave.cpp)
target_link_libraries(test_wave PRIVATE GTest::gtest_main mfem-gbs)
target_include_directories(test_wave PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_definitions(test_wave PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# alfven app
add_executable(test_alfven unit/test_alfven.cpp)
target_link_libraries(test_alfven PRIVATE GTest::gtest_main mfem-gbs)
target_include_directories(test_alfven PRIVATE ${CMAKE_SOURCE_DIR}/src)
target_compile_definitions(test_alfven PRIVATE DATA_DIR="${CMAKE_SOURCE_DIR}/data")

# Discover and register tests
include(GoogleTest)

gtest_discover_tests(test_wave   WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )
gtest_discover_tests(test_alfven WORKING_DIRECTORY ${CMAKE_BINARY_DIR} )

