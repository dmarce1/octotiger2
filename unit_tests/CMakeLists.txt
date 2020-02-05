# Copyright (c) 2015 Thomas Heller
# Copyright (c) 2015 Dominic Marcello
# Copyright (c) 2018-2019 Parsa Amini
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)


add_hpx_executable(
  super_array_test
  DEPENDENCIES
    Silo::silo Boost::boost 
  SOURCES
    super_array.cpp
 )


target_include_directories(super_array_test PUBLIC ${PROJECT_SOURCE_DIR})

