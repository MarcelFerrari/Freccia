cmake_minimum_required(VERSION 3.10)
project(Freccia)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wno-enum-compare -fopenmp -march=native -std=c++17 -DEIGEN_NO_DEBUG -DEIGEN_USE_MKL_ALL -DMKL_ILP64 -m64 -DFRECCIA_USE_MKL")

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories($ENV{MKLROOT}/include)

# Linker flags
link_directories($ENV{MKLROOT}/lib/intel64)

# Source files
file(GLOB_RECURSE SOURCES "src/*.cpp")

# Enable position independent code - required for linking with Python
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Create the Freccia static library
add_library(Freccia STATIC ${SOURCES})

# Set the output directory for the static library
set_target_properties(Freccia PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/lib")

# Linking MKL and other required libraries
target_link_libraries(Freccia mkl_intel_ilp64 mkl_tbb_thread mkl_core pthread m dl)

# Installation of library
install(TARGETS Freccia DESTINATION lib)
