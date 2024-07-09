cmake_minimum_required(VERSION 3.15)
project(affine)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(GTest REQUIRED)
find_library(CAPD_LIBRARY capd PATHS $ENV{HOME}/.local/lib)

enable_testing()

add_compile_options(-Wall -Wextra -Wpedantic)

add_executable(affine_tests)

target_sources(
    affine_tests
    PUBLIC
    tests/constructor_test.cpp
    tests/element_test.cpp
    tests/spanning_space_test.cpp
    tests/intersection_test.cpp
)

target_include_directories(
    affine_tests
    PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    $ENV{HOME}/.local/include
)

target_link_libraries(
    affine_tests
    PUBLIC
    ${CAPD_LIBRARY}
    GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(affine_tests)