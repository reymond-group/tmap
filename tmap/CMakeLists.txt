cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(MSVC)
    set(CMAKE_GENERATOR_TOOLSET "ClangCL")
endif(MSVC)

project(_tmap)

if(MSVC)
    # set(CMAKE_CXX_FLAGS_DEBUG "/O2 /permissive /std:c++17")
    # set(CMAKE_CXX_FLAGS_RELEASE "/O2 /permissive /std:c++17 /wd4068 /wd4305 /wd4267 /wd4244")
else()
    set(CMAKE_CXX_FLAGS_DEBUG "-O3 -funroll-loops -fopenmp")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -fopenmp")
endif(MSVC)


SET(OGDF_INCLUDE "/include")
SET(OGDF_LIB "/lib")

include_directories(${OGDF_INCLUDE})
link_directories(${OGDF_LIB})

# find_package(OpenMP)
# if (OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
# endif()

find_package(OpenMP REQUIRED)

set(THREADS_PREFER_PTHREAD_FLAG ON)

find_package(Threads REQUIRED)
# include_directories(${pybind11_INCLUDE_DIRS})

SET(SOURCE_DIR "src/_tmap")
# Tell cmake that headers are in alse in source_dir
include_directories(${SOURCE_DIR})
SET(SOURCES "${SOURCE_DIR}/layout.cc" "${SOURCE_DIR}/lshforest.cc" "${SOURCE_DIR}/minhash.cc" "${SOURCE_DIR}/analyse.cc")

# Generate python module
# find_package(pybind11 CONFIG REQUIRED)
add_subdirectory(pybind11)

pybind11_add_module(_tmap ${SOURCES} "${SOURCE_DIR}/bindings.cc")
if(MSVC)
target_link_libraries(_tmap PRIVATE OGDF.lib COIN.lib)
else()
target_link_libraries(_tmap PRIVATE OGDF COIN)
endif(MSVC)

target_link_libraries(_tmap PRIVATE Threads::Threads)
target_link_libraries(_tmap PRIVATE ${OpenMP_CXX_FLAGS})
