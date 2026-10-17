cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(channel_with_leaflet_ref C CXX Fortran)

# Point OOMPH_INSTALL at an oomph-lib installation, e.g.
#   cmake -G Ninja -B build -DOOMPH_INSTALL=$HOME/oomph-lib/install
find_package(oomphlib CONFIG REQUIRED PATHS "${OOMPH_INSTALL}")

oomph_add_executable(
  NAME channel_with_leaflet_ref
  SOURCES channel_with_leaflet_ref.cc
  LIBRARIES oomph::multi_physics oomph::beam oomph::navier_stokes
            oomph::meshes oomph::generic)
