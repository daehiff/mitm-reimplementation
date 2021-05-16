cmake_minimum_required(VERSION 3.14)
project(mitm_algorithm)

set(CMAKE_CXX_STANDARD 14)

find_package(xtl REQUIRED)
find_package(xtensor REQUIRED)
find_package(xtensor-blas REQUIRED)
find_package(xtensor-python REQUIRED)

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


set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -Wall -Wextra -Werror")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -Wall -Wextra -Werror")



add_executable(mitm_algorithm
        src/mitm.cpp
        src/circuit.cpp
        src/circuit.h
        src/gate.cpp
        src/gate.h
        src/utils.h
        src/dispatch_queue.h)

target_link_libraries(mitm_algorithm xtensor xtensor-blas xtensor-python)



