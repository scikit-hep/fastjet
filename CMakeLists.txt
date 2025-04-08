cmake_minimum_required(VERSION 3.11)

if(NOT SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")  # provided by scikit-build-core
endif()
string(REPLACE "." ";" VERSION_SPLIT ${SKBUILD_PROJECT_VERSION})
list(GET VERSION_SPLIT 0 SPLIT_VERSION_MAJOR)
list(GET VERSION_SPLIT 1 SPLIT_VERSION_MINOR)

project(fastjet VERSION ${SPLIT_VERSION_MAJOR}.${SPLIT_VERSION_MINOR} LANGUAGES CXX)

set(PKG_INSTALL "fastjet") ## python package name

find_package(pybind11 REQUIRED)

add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-core)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-contrib)

pybind11_add_module(_ext MODULE src/_ext.cpp)
target_link_libraries(_ext PRIVATE fastjet::fastjet)
target_link_libraries(_ext PRIVATE EnergyCorrelator LundPlane Nsubjettiness RecursiveTools)
target_include_directories(_ext PRIVATE ${CMAKE_BINARY_DIR}/include)
#set_target_properties(_ext PROPERTIES BUILD_WITH_INSTALL_RPATH ON)

install(TARGETS _ext LIBRARY DESTINATION ${PKG_INSTALL})
