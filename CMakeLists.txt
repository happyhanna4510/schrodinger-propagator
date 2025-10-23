cmake_minimum_required(VERSION 3.15)
project(schrodinger_propagator LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(CTest)

set(SOURCES
    src/cli.cpp
    src/grid.cpp
    src/core/io_utils.cpp
    src/core/math_utils.cpp
    src/core/spectral.cpp
    src/evolve/evolve_factory.cpp
    src/evolve/evolve_rk4.cpp
    src/evolve/evolve_taylor.cpp
    src/hamiltonian.cpp
    src/initial.cpp
    src/io.cpp
    src/main.cpp
    src/morse_analytic.cpp
    src/morse_potential.cpp
    src/morse_static.cpp
    src/paths.cpp
    src/snapshots.cpp
    src/taylor.cpp)

add_executable(morse ${SOURCES})

target_include_directories(morse PRIVATE ${CMAKE_SOURCE_DIR}/include)

target_compile_features(morse PRIVATE cxx_std_17)

if (MSVC)
    target_compile_options(morse PRIVATE /W4 /permissive-)
else()
    target_compile_options(morse PRIVATE -Wall -Wextra -Wpedantic)
endif()

find_package(Eigen3 3.3 QUIET CONFIG)
if (NOT Eigen3_FOUND)
    find_package(Eigen3 3.3 QUIET)
endif()
if (NOT Eigen3_FOUND)
    set(EIGEN_FALLBACK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/eigen-5.0.0")
    if (EXISTS "${EIGEN_FALLBACK_DIR}/Eigen/Dense")
        add_library(Eigen3::Eigen INTERFACE IMPORTED)
        target_include_directories(Eigen3::Eigen INTERFACE "${EIGEN_FALLBACK_DIR}")
    else()
        message(FATAL_ERROR "Eigen3 not found. Install it or provide third_party/eigen-5.0.0")
    endif()
endif()

target_link_libraries(morse PRIVATE Eigen3::Eigen)

if (BUILD_TESTING)
    set(TEST_OUTPUT_DIR "${CMAKE_BINARY_DIR}/test-artifacts")
    add_test(NAME morse_smoke
        COMMAND morse
            --gamma 5
            --N 101
            --xmax 10
            --evolve taylor
            --dt 1e-4
            --tmax 1e-4
            --K 4
            --log 1
            --no_wide
            --quiet
            --outdir ${TEST_OUTPUT_DIR})
    set_tests_properties(morse_smoke PROPERTIES WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
endif()
