cmake_minimum_required(VERSION 3.15)
project(MyProject LANGUAGES C Fortran)

find_package(Python COMPONENTS Interpreter Development REQUIRED)
#find_package(scikit-build REQUIRED)
find_package(Cython REQUIRED)

# Install Python scripts
install(DIRECTORY src/ase2sprkkr DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages)

add_subdirectory(src/ase2sprkkr/bindings/xband)
