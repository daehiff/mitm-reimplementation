cmake_minimum_required(VERSION 3.14)
project(mitm_algorithm)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)


find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)

if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS /MANIFEST:NO)
endif ()

add_definitions(-DHAVE_CBLAS=1)

if (WIN32)
    find_package(OpenBLAS REQUIRED)
    set(BLAS_LIBRARIES ${CMAKE_INSTALL_PREFIX}${OpenBLAS_LIBRARIES})
else ()
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
endif ()

message(STATUS "BLAS VENDOR:    " ${BLA_VENDOR})
message(STATUS "BLAS LIBRARIES: " ${BLAS_LIBRARIES})


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

add_executable(mitm_algorithm main.cpp Circuit.cpp Circuit.h Gate.cpp Gate.h utils.h)


target_link_libraries(mitm_algorithm xtensor xtensor-blas ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

target_compile_options(mitm_algorithm PRIVATE -march=native)

