# Copyright (c) 2015 Thomas Heller
# Copyright (c) 2015 Dominic Marcello
# Copyright (c) 2018-2019 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

cmake_minimum_required(VERSION 3.12.4)
cmake_policy(VERSION 3.12.4)

# Libraries linked via full path no longer produce linker search paths.
cmake_policy(SET CMP0003 NEW)
cmake_policy(SET CMP0074 NEW)

project(octotiger2 CXX)

# Search path for CMake modules to be loaded by include() and find_package()
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

find_package(HPX REQUIRED NO_CMAKE_PACKAGE_REGISTRY)
find_package(Silo REQUIRED)
find_package(Boost REQUIRED)

set(octotiger2_source_files
	src/conserved.cpp
	src/main.cpp
	src/options.cpp
	src/primitive.cpp
    src/riemann.cpp
	src/tree.cpp
)

set(octotiger2_header_files
    octotiger/dim.hpp
    octotiger/conserved.hpp
    octotiger/fixed_real.hpp
    octotiger/math.hpp
	octotiger/options.hpp
	octotiger/primitive.hpp
    octotiger/real.hpp
    octotiger/riemann.hpp
	octotiger/super_array.hpp
	octotiger/tree.hpp
    octotiger/vect.hpp
    octotiger/volume.hpp
)

add_hpx_executable(
  octotiger2
  DEPENDENCIES
    Silo::silo Boost::boost 
  SOURCES
    ${octotiger2_source_files}
  HEADERS
    ${octotiger2_header_files}
 )


target_include_directories(octotiger2 PUBLIC ${PROJECT_SOURCE_DIR})
add_subdirectory(unit_tests)