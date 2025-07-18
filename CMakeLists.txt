cmake_minimum_required(VERSION 3.15...4.0)

if(NOT SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")  # provided by scikit-build-core
endif()
string(REPLACE "." ";" VERSION_SPLIT ${SKBUILD_PROJECT_VERSION})
list(GET VERSION_SPLIT 0 SPLIT_VERSION_MAJOR)
list(GET VERSION_SPLIT 1 SPLIT_VERSION_MINOR)

project(fastjet VERSION ${SPLIT_VERSION_MAJOR}.${SPLIT_VERSION_MINOR} LANGUAGES CXX)

option(SKHEPFJ_USE_INSTALLED_FASTJET         "Use an existing installed version of fastjet"                         OFF)
option(SKHEPFJ_USE_INSTALLED_FASTJET_CONTRIB "Use an existing installed version of fastjet-contrib"                 OFF)
option(SKHEPFJ_PATCH_FJCORE                  "Patch fastjet-contrib with CMAKE_INSTALL_LIBDIR use"                  OFF)
option(SKHEPFJ_PATCH_FJCONTRIB               "Patch fastjet-contrib with cmake build"                               ON)
option(SKHEPFJ_PASSTHRU_FASTJET_SWIG         "Create a _swig.py that just imports <fastjet-lib python module>"      OFF)
option(SKHEPFJ_COPY_GMP_DLL                  "On Windows, force a copy of the amp runtime .dll to the wheel area"   OFF)

set(SKHEPFJ_FASTJET_LIB_PYMODULE_NAME "fastjet_cxx" CACHE STRING "The name of the installed fastjet python package.")

set(PKG_INSTALL "fastjet") ## python package name

include(FetchContent)

find_package(pybind11 REQUIRED)

if (SKHEPFJ_USE_INSTALLED_FASTJET)
  find_package(fastjet REQUIRED)
else()
  add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-core)
endif()

if (SKHEPFJ_USE_INSTALLED_FASTJET_CONTRIB)
  find_package(fastjetcontrib REQUIRED)
else()
  add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-contrib)
endif()

if (SKHEPFJ_PASSTHRU_FASTJET_SWIG)
  file(WRITE ${CMAKE_INSTALL_PREFIX}/_swig.py "from ${SKHEPFJ_FASTJET_LIB_PYMODULE_NAME} import *\n")
endif()

pybind11_add_module(_ext MODULE src/_ext.cpp)
target_link_libraries(_ext
                      PRIVATE
                        fastjet::fastjet
                        fastjet::contrib::fastjetcontribfragile
                     )
if (MSVC)
  target_compile_definitions(_ext PRIVATE _USE_MATH_DEFINES)
endif()
target_include_directories(_ext PRIVATE ${CMAKE_BINARY_DIR}/include)
#set_target_properties(_ext PROPERTIES BUILD_WITH_INSTALL_RPATH TRUE)
#set_target_properties(_ext PROPERTIES INSTALL_RPATH_USE_LINK_PATH TRUE)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set_target_properties(_ext PROPERTIES INSTALL_RPATH "@loader_path/${CMAKE_INSTALL_LIBDIR};@loader_path/../../../../${CMAKE_INSTALL_LIBDIR}")
elseif(NOT MSVC)
  set_target_properties(_ext PROPERTIES INSTALL_RPATH "$ORIGIN/${CMAKE_INSTALL_LIBDIR}:$ORIGIN/../../../../${CMAKE_INSTALL_LIBDIR}")
endif()

if (SKHEPFJ_COPY_GMP_DLL AND MSVC AND FASTJET_ENABLE_CGAL)
  message("-- Copying ${GMP_DLL_TO_COPY} to project area!")
  install(FILES ${GMP_DLL_TO_COPY} DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
install(TARGETS _ext LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
