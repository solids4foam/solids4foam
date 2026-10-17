cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
project(collapsible_channel_ref C CXX Fortran)

# Point OOMPH_INSTALL at an oomph-lib installation, e.g.
#   cmake -G Ninja -B build -DOOMPH_INSTALL=$HOME/oomph-lib/install
find_package(oomphlib CONFIG REQUIRED PATHS "${OOMPH_INSTALL}")

oomph_add_executable(
  NAME collapsible_channel_ref_CR
  SOURCES collapsible_channel_ref.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic)

oomph_add_executable(
  NAME collapsible_channel_ref_TH
  SOURCES collapsible_channel_ref.cc
  LIBRARIES oomph::beam oomph::navier_stokes oomph::meshes oomph::generic
  CXX_DEFINITIONS TAYLOR_HOOD)
