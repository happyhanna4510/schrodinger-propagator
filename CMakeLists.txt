cmake_minimum_required(VERSION 3.15)
project(morse_demo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Где лежат заголовки Eigen и твои include/
# ВАЖНО: путь должен содержать папку Eigen/, т.е. third_party/eigen-5.0.0/Eigen/Core и т.п.
set(EIGEN_ROOT ${CMAKE_SOURCE_DIR}/third_party/eigen-5.0.0)

set(SRC
    src/main.cpp
    src/grid.cpp              # <<< добавь
    src/hamiltonian.cpp
    src/morse_potential.cpp
    src/morse_analytic.cpp    # <<< добавь
    src/io.cpp
    src/taylor.cpp
    src/cli.cpp
    src/paths.cpp
    src/initial.cpp
    src/snapshots.cpp
    src/morse_static.cpp
)

add_executable(morse ${SRC})

# Подключаем заголовки к цели (вместо глобального include_directories)
target_include_directories(morse PRIVATE
    ${EIGEN_ROOT}
    ${CMAKE_SOURCE_DIR}/include
)

# НИЧЕГО НЕ ЛИНКУЕМ С Eigen — это header-only
# target_link_libraries(morse PRIVATE Eigen3::Eigen)  # <-- УДАЛИТЬ

# На всякий случай
if (MSVC)
  target_compile_options(morse PRIVATE /std:c++17)
endif()
