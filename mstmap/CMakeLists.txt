cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 14)
# set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O3 -funroll-loops")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops")

project(tmap)


if (DEFINED ENV{CONDA_PREFIX})
    SET(OGDF_INCLUDE "$ENV{CONDA_PREFIX}/include")
    SET(OGDF_LIB "$ENV{CONDA_PREFIX}/lib")

    include_directories(${OGDF_INCLUDE})
    link_directories(${OGDF_LIB})
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

SET(SOURCE_DIR "tmap")
# Tell cmake that headers are in alse in source_dir
include_directories(${SOURCE_DIR})
SET(SOURCES "${SOURCE_DIR}/layout.cc" "${SOURCE_DIR}/lshforest.cc" "${SOURCE_DIR}/minhash.cc")

# Generate python module
add_subdirectory(lib/pybind11)
pybind11_add_module(tmap ${SOURCES} "${SOURCE_DIR}/bindings.cc")
target_link_libraries(tmap PRIVATE OGDF COIN)
target_link_libraries(tmap PRIVATE Threads::Threads)