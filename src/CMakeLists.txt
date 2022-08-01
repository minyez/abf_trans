# Specify the executable target
add_executable(abf_trans_exe main.cpp)
# The default name of the executable will be exactly the variable name
# To change the varaiable, we use `set_target_properties`
set_target_properties(abf_trans_exe PROPERTIES
  OUTPUT_NAME "abf_trans.exe")

# Now specify where the executable will be generated.
# By default, if the compilation succeeds, you will have the executable at the "src" subdirectory in the build directory
# To put it in the build root directory
set_target_properties(abf_trans_exe
  PROPERTIES
  ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}
  RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})

# the source files to use
list(APPEND _sources
  abf.cpp
  cell.cpp
  # constants.cpp
  io.cpp
  kmesh.cpp
  rotate.cpp
  spglib_utils.cpp
  sph.cpp)
# set include directories, required when external C/C++ codes are used in the project
target_include_directories(abf_trans_exe
  PUBLIC
    ${SPGLIB_PATH}/include
  )

target_sources(abf_trans_exe PRIVATE ${_sources})
target_link_libraries(abf_trans_exe ${LIBS})
