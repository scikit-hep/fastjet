cmake_minimum_required(VERSION 3.15...4.0)

if(NOT SKBUILD_PROJECT_VERSION)
  set(SKBUILD_PROJECT_VERSION "0.0.0")  # provided by scikit-build-core
endif()
string(REPLACE "." ";" VERSION_SPLIT ${SKBUILD_PROJECT_VERSION})
list(GET VERSION_SPLIT 0 SPLIT_VERSION_MAJOR)
list(GET VERSION_SPLIT 1 SPLIT_VERSION_MINOR)

project(fastjet VERSION ${SPLIT_VERSION_MAJOR}.${SPLIT_VERSION_MINOR} LANGUAGES CXX)

option(SKHEPFJ_USE_INSTALLED_FASTJET          "Use an existing installed version of fastjet"         OFF)
option(SKHEPFJ_USE_INSTALLED_FASTJET_CONTRIB  "Use an existing installed version of fastjet-contrib" OFF)
option(SKHEPFJ_PATCH_FJCORE                   "Patch fastjet-contrib with CMAKE_INSTALL_LIBDIR use"  ON)
option(SKHEPFJ_PATCH_FJCONTRIB                "Patch fastjet-contrib with cmake build"               ON)

set(PKG_INSTALL "fastjet") ## python package name

find_package(pybind11 REQUIRED)

if (SKHEPFJ_USE_INSTALLED_FASTJET)
  find_package(fastjet REQUIRED)
else()
  if (SKHEPFJ_PATCH_FJCORE)
    execute_process(COMMAND
                      git apply --whitespace=fix ../../0001-use-CMAKE_INSTALL_LIBDIR-consistently-in-pyinterface.patch
                    WORKING_DIRECTORY
                      ${CMAKE_CURRENT_SOURCE_DIR}/extern/fastjet-core
                   )
  endif()
  add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-core)
endif()

if (SKHEPFJ_USE_INSTALLED_FASTJET_CONTRIB)
  find_package(fastjetcontrib REQUIRED)
else()
  if (SKHEPFJ_PATCH_FJCONTRIB)
    execute_process(COMMAND
                      git apply --whitespace=fix ../../0001-fjcontrib-top-level-only-cmake-scripts.patch
                    WORKING_DIRECTORY
                      ${CMAKE_CURRENT_SOURCE_DIR}/extern/fastjet-contrib
                   )
  endif()
  add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-contrib)
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
  set_target_properties(_ext PROPERTIES INSTALL_RPATH "@loader_path/${CMAKE_INSTALL_LIBDIR}")
elseif(NOT MSVC)
  set_target_properties(_ext PROPERTIES INSTALL_RPATH "$ORIGIN/${CMAKE_INSTALL_LIBDIR}")
endif()

if (MSVC AND FASTJET_ENABLE_CGAL)
  message("-- Copying ${GMP_DLL_TO_COPY} to project area!")
  install(FILES ${GMP_DLL_TO_COPY} DESTINATION ${CMAKE_INSTALL_PREFIX})
endif()
install(TARGETS _ext LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX})
