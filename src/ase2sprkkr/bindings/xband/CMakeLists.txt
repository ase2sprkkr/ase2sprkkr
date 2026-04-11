# Base Fortran flags
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy -finit-integer=0")

# Debug-specific options
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy -O0 -finit-local-zero -g -fcheck=all")
  add_compile_options(-g)
endif()

find_package(Python COMPONENTS Interpreter Development NumPy REQUIRED)
include_directories(${Python_NumPy_INCLUDE_DIRS})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB SYMMETRY_SOURCES "symmetry/*.f")
add_library(symmetry_lib STATIC ${SYMMETRY_SOURCES})
target_compile_definitions(symmetry_lib PUBLIC MYLIB_STATIC_DEFINE)
set_target_properties(symmetry_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)


add_custom_command(
  OUTPUT symmetry.c
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/symmetry.c from ${CMAKE_CURRENT_SOURCE_DIR}/symmetry.pyx"
  COMMAND Python::Interpreter -m cython
          "${CMAKE_CURRENT_SOURCE_DIR}/symmetry.pyx" --output-file symmetry.c
  DEPENDS symmetry.pyx
  VERBATIM)

python_add_library(symmetry MODULE symmetry.c WITH_SOABI)
target_link_libraries(symmetry PRIVATE symmetry_lib)

install(TARGETS symmetry
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/ase2sprkkr/bindings/xband)


set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finit-integer=0")
# Include both fixed-form and free-form Fortran sources
file(GLOB SPHERES_F_SOURCES "spheres/*.f90")
file(GLOB SPHERES_F90_SOURCES "spheres/*.f90")
add_library(spheres_lib STATIC ${SPHERES_F90_SOURCES})
target_compile_definitions(spheres_lib PUBLIC MYLIB_STATIC_DEFINE)
set_target_properties(spheres_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Ensure free-form Fortran is recognized where needed
foreach(src ${SPHERES_F90_SOURCES})
  set_source_files_properties(${src} PROPERTIES Fortran_FORMAT FREE)
endforeach()

add_custom_command(
  OUTPUT spheres.c
  COMMENT
    "Making ${CMAKE_CURRENT_BINARY_DIR}/spheres.c from ${CMAKE_CURRENT_SOURCE_DIR}/spheres.pyx"
  COMMAND Python::Interpreter -m cython
          "${CMAKE_CURRENT_SOURCE_DIR}/spheres.pyx" --output-file spheres.c
  DEPENDS spheres.pyx
  VERBATIM)

python_add_library(spheres MODULE spheres.c WITH_SOABI)
target_link_libraries(spheres PRIVATE spheres_lib)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  target_compile_options(spheres PRIVATE -g)
  target_link_options(spheres PRIVATE -g)
  set(CMAKE_STRIP FALSE)
  set(CMAKE_INSTALL_DO_STRIP FALSE)
endif()

install(TARGETS spheres
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/ase2sprkkr/bindings/xband)

