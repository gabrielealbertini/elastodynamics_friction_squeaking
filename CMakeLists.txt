#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
# @section DESCRIPTION
#
#===============================================================================

cmake_minimum_required(VERSION 2.6)
enable_language(CXX)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

project(akantu_simulations)

find_package(Akantu REQUIRED)
find_package(GMSH REQUIRED)
find_package(uguca)

find_package(Boost REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

if(AKANTU_CORE_CXX11)
  message(TOTO)
  add_definition(-std=c++11)
endif()

include(cmake/AkantuSimulationsMacros.cmake)

include_directories(${AKANTU_INCLUDE_DIRS}
		    ${AKANTU_BOOST_INCLUDE_DIR})


if(CMAKE_CXX_FLAGS STREQUAL "")
  set(CMAKE_CXX_FLAGS "-Wall" CACHE STRING "Flags used by the compiler during all build types." FORCE)
endif()

# Debug
#option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
#mark_as_advanced(AKANTU_DEBUG)

add_definitions(-DAKANTU_NDEBUG)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -DAKANTU_NDEBUG"
  CACHE STRING "Flags used by the compiler during release builds" FORCE)

set(CMAKE_CXX_FLAGS_PROFILING "-g -pg" CACHE STRING "Flags used by the compiler during profiling builds")
set(CMAKE_EXE_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
set(CMAKE_SHARED_LINKER_FLAGS_PROFILING "-pg" CACHE STRING "Flags used by the linker during profiling builds")
mark_as_advanced(CMAKE_CXX_FLAGS_PROFILING)
mark_as_advanced(CMAKE_EXE_LINKER_FLAGS_PROFILING)
mark_as_advanced(CMAKE_SHARED_LINKER_FLAGS_PROFILING)

#===============================================================================
# List of simulations
#===============================================================================

add_simulation(bimat_friction "bimaterial friction - Adam's instability")

#add_simulation(crack_path "crack path at interfaces")

#add_simulation(friction_fracture_comparison "Comparison of friction and fracture for Ilya's simulations")

#add_simulation(hybrid_method "Hybrid FEM SBI method for simulating infinite boundary")

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/local_simulations.cmake")
   include("local_simulations.cmake")
endif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/local_simulations.cmake")

