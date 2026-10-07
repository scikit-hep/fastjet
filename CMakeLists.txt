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

# Before add_subdirectory so the Python targets (incl. SABIModule in abi3 builds) are visible to the SWIG module.
find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)

if (SKHEPFJ_USE_INSTALLED_FASTJET)
  find_package(fastjet REQUIRED)
else()
  add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-core)
endif()

if (TARGET fastjet_swig)
  # PyUnicode_AsUTF8 is outside the 3.12 limited API; a patched copy of the header shadows the submodule's.
  set(_fjpyext ${CMAKE_SOURCE_DIR}/extern/fastjet-core/pyinterface/FastJetPythonExtensions.hh)
  file(READ ${_fjpyext} _fjpyext_src)
  string(REPLACE "PyUnicode_AsUTF8(py_str)" "PyUnicode_AsUTF8AndSize(py_str, NULL)" _fjpyext_src "${_fjpyext_src}")
  file(CONFIGURE OUTPUT ${CMAKE_BINARY_DIR}/swig-shim/FastJetPythonExtensions.hh CONTENT "${_fjpyext_src}" @ONLY)
  set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_fjpyext})
  target_include_directories(fastjet_swig BEFORE PRIVATE ${CMAKE_BINARY_DIR}/swig-shim)

  if (SKBUILD_SABI_COMPONENT)
    # SWIG < 4.3 output calls malloc/free without <cstdlib>, which Python.h only includes outside the limited API.
    if (SWIG_VERSION VERSION_LESS 4.3)
      message(FATAL_ERROR "SWIG ${SWIG_VERSION} cannot build the limited-API (abi3) module: install SWIG >= 4.3, "
                          "or build a version-specific wheel with -Cwheel.py-api=")
    endif()
    # 3.12 is the floor nanobind's STABLE_ABI imposes on _ext, so the wheel tag holds for both modules.
    target_compile_definitions(fastjet_swig PRIVATE Py_LIMITED_API=0x030C0000)
    if (MSVC)
      set_target_properties(fastjet_swig PROPERTIES SUFFIX ".pyd")
    else()
      set_target_properties(fastjet_swig PROPERTIES SUFFIX ".abi3${CMAKE_SHARED_MODULE_SUFFIX}")
    endif()
    foreach(_prop LINK_LIBRARIES INTERFACE_LINK_LIBRARIES)
      get_target_property(_libs fastjet_swig ${_prop})
      list(TRANSFORM _libs REPLACE "^Python::Module$" "Python::SABIModule")
      set_target_properties(fastjet_swig PROPERTIES ${_prop} "${_libs}")
    endforeach()
  endif()
endif()

if (SKHEPFJ_USE_INSTALLED_FASTJET_CONTRIB)
  find_package(fastjetcontrib REQUIRED)
else()
  add_subdirectory(${CMAKE_SOURCE_DIR}/extern/fastjet-contrib)
endif()

if (SKHEPFJ_PASSTHRU_FASTJET_SWIG)
  file(WRITE ${CMAKE_INSTALL_PREFIX}/_swig.py "from ${SKHEPFJ_FASTJET_LIB_PYMODULE_NAME} import *\n")
endif()

# NB_DOMAIN keeps the catch-all exception translator private to fastjet; MUSL_DYNAMIC_LIBCPP shares libstdc++ with libfastjet.
nanobind_add_module(_ext STABLE_ABI NB_DOMAIN fastjet MUSL_DYNAMIC_LIBCPP src/_ext.cpp)
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
