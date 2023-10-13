cmake_minimum_required(VERSION 3.10)
project(FrecciaTests)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-enum-compare -fopenmp -march=native -std=c++17 -DEIGEN_NO_DEBUG -DEIGEN_USE_MKL_ALL -DMKL_ILP64 -m64")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/../include/)
include_directories($ENV{MKLROOT}/include)

# Linker flags
link_directories($ENV{MKLROOT}/lib/intel64)
link_directories(${CMAKE_SOURCE_DIR}/../lib)

# For every source file, create an executable and link it.
file(GLOB TEST_SOURCES "src/*.cpp")
foreach(testsourcefile ${TEST_SOURCES})
    # Extract the filename without extension
    get_filename_component(testname ${testsourcefile} NAME_WE)
    
    # Add an executable for the source file
    add_executable(${testname} ${testsourcefile})

    # Link the executable with Freccia and MKL libraries
    target_link_libraries(${testname} 
        Freccia
        mkl_intel_ilp64
        mkl_tbb_thread 
        mkl_core
        pthread
        m
        dl
    )
endforeach(testsourcefile ${TEST_SOURCES})
