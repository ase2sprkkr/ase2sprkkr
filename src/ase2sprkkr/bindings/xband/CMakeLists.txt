set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -std=legacy")
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



file(GLOB SPHERES_SOURCES "spheres/*.f")
add_library(spheres_lib STATIC ${SPHERES_SOURCES})
target_compile_definitions(spheres_lib PUBLIC MYLIB_STATIC_DEFINE)
set_target_properties(spheres_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

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

install(TARGETS spheres
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/ase2sprkkr/bindings/xband)


