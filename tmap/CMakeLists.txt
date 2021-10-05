cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)

project(tmap)

if(MSVC)
    set(CMAKE_CXX_FLAGS_DEBUG "/openmp /O2")
    set(CMAKE_CXX_FLAGS_RELEASE "/openmp /O2")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-O3 -funroll-loops -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -fopenmp")
endif(MSVC)


if (DEFINED ENV{CONDA_PREFIX})
    SET(OGDF_INCLUDE "$ENV{CONDA_PREFIX}/include")
    SET(OGDF_LIB "$ENV{CONDA_PREFIX}/lib")

    include_directories(${OGDF_INCLUDE})
    link_directories(${OGDF_LIB})
endif()

# find_package(OpenMP)
# if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

find_package(OpenMP REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

SET(SOURCE_DIR "tmap")
# Tell cmake that headers are in alse in source_dir
include_directories(${SOURCE_DIR})
SET(SOURCES "${SOURCE_DIR}/layout.cc" "${SOURCE_DIR}/lshforest.cc" "${SOURCE_DIR}/minhash.cc" "${SOURCE_DIR}/analyse.cc")

# Generate python module
add_subdirectory(lib/pybind11)

pybind11_add_module(tmap ${SOURCES} "${SOURCE_DIR}/bindings.cc")
if(MSVC)
target_link_libraries(tmap PRIVATE OGDF.lib COIN.lib)
else()
target_link_libraries(tmap PRIVATE OGDF COIN)
endif(MSVC)

target_link_libraries(tmap PRIVATE Threads::Threads)
target_link_libraries(tmap PRIVATE ${OpenMP_CXX_FLAGS})
