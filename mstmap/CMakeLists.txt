cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(mstmap)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

SET(SOURCE_DIR "src/mstmap")
# Tell cmake that headers are in alse in source_dir
include_directories(${SOURCE_DIR})
SET(SOURCES "${SOURCE_DIR}/layout.cc" "${SOURCE_DIR}/lshforest.cc")

# Generate python module
add_subdirectory(lib/pybind11)
pybind11_add_module(mstmap ${SOURCES} "${SOURCE_DIR}/bindings.cc")
target_link_libraries(mstmap PRIVATE OGDF COIN)
target_link_libraries(mstmap PRIVATE Threads::Threads)