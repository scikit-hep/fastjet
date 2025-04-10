cmake_minimum_required(VERSION 3.15...4.0)

if(NOT SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")  # provided by scikit-build-core
endif()
string(REPLACE "." ";" VERSION_SPLIT ${SKBUILD_PROJECT_VERSION})
list(GET VERSION_SPLIT 0 SPLIT_VERSION_MAJOR)
list(GET VERSION_SPLIT 1 SPLIT_VERSION_MINOR)

project(fastjet VERSION ${SPLIT_VERSION_MAJOR}.${SPLIT_VERSION_MINOR} LANGUAGES CXX)

set(PKG_INSTALL "fastjet") ## python package name

include(FetchContent)

find_package(pybind11 REQUIRED)

find_package(CGAL)
if (NOT CGAL_FOUND)
  FetchContent_Declare(CGAL
    GIT_REPOSITORY https://github.com/CGAL/cgal.git
    GIT_TAG  v6.0.1
  )
  FetchContent_MakeAvailable(CGAL)
endif()

#FetchContent_Declare(fastjet
#    SOURCE_DIR ${CMAKE_SOURCE_DIR}/extern/fastjet-core
#    CMAKE_ARGS "-DCMAKE_BUILD_TYPE=\"Release\" -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/${PKG_INSTALL} -DFASTJET_CUSTOM_PYTHON_INSTALL=lib/${PKG_INSTALL} -DFASTJET_PYTHON_PACKAGE_NAME=_swig -DFASTJET_ENABLE_CGAL=ON"
#  )
#FetchContent_MakeAvailable(fastjet)

add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-core)
add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-contrib)

pybind11_add_module(_ext MODULE src/_ext.cpp)
target_link_libraries(_ext PRIVATE fastjet::fastjet)
target_link_libraries(_ext PRIVATE EnergyCorrelator LundPlane Nsubjettiness RecursiveTools)
target_include_directories(_ext PRIVATE ${CMAKE_BINARY_DIR}/include)
#set_target_properties(_ext PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
#set_target_properties(_ext PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(_ext PROPERTIES INSTALL_RPATH "@loader_path/lib")
else()
  set_target_properties(_ext PROPERTIES INSTALL_RPATH "$ORIGIN/lib")
endif()

install(TARGETS _ext LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
