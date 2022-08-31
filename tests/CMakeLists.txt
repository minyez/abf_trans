add_executable(test_sph test_sph.cpp)
target_sources(test_sph PRIVATE ../src/sph.cpp)

add_executable(test_linalg test_linalg.cpp)
target_link_libraries(test_linalg "${LAPACK_LIBRARIES}")

add_executable(test_matrix test_matrix.cpp)
target_link_libraries(test_matrix "${LAPACK_LIBRARIES}")

add_executable(test_spglib_utils ./test_spglib_utils.cpp)
target_sources(test_spglib_utils PRIVATE ../src/spglib_utils.cpp)
target_link_libraries(test_spglib_utils "${SPGLIB_LIBRARIES} ${LAPACK_LIBRARIES}")

add_executable(test_abf ./test_abf.cpp)
target_sources(test_abf PRIVATE ../src/abf.cpp)

add_executable(test_rotate test_rotate.cpp)
target_sources(test_rotate PRIVATE ../src/rotate.cpp)
target_link_libraries(test_rotate "${LAPACK_LIBRARIES}")

add_executable(test_mathtools ./test_mathtools.cpp)

add_executable(test_kmesh test_kmesh.cpp)
target_sources(test_kmesh PRIVATE ../src/kmesh.cpp ../src/cell.cpp ../src/spglib_utils.cpp)
target_link_libraries(test_kmesh "${SPGLIB_LIBRARIES} ${LAPACK_LIBRARIES}")

add_executable(test_vec test_vec.cpp)
target_link_libraries(test_vec "${LAPACK_LIBRARIES}")

add_executable(test_cell test_cell.cpp)
target_sources(test_cell PRIVATE ../src/cell.cpp)
target_link_libraries(test_cell "${LAPACK_LIBRARIES}")

add_executable(test_trans test_trans.cpp)
target_sources(test_trans PRIVATE ../src/trans.cpp ../src/logger.cpp ../src/kmesh.cpp ../src/abf.cpp ../src/cell.cpp ../src/rotate.cpp ../src/spglib_utils.cpp)
target_link_libraries(test_trans "${SPGLIB_LIBRARIES} ${LAPACK_LIBRARIES}")

add_executable(test_io test_io.cpp)
target_sources(test_io PRIVATE ../src/io.cpp ../src/abf.cpp)
target_link_libraries(test_io "${LAPACK_LIBRARIES}")
# Copy sample files for input test
configure_file(sample_complexmat.csc ${CMAKE_CURRENT_BINARY_DIR}/sample_complexmat.csc COPYONLY)


add_test(
  NAME test_sph
  COMMAND $<TARGET_FILE:test_sph>
  )

add_test(
  NAME test_linalg
  COMMAND $<TARGET_FILE:test_linalg>
  )

add_test(
  NAME test_matrix
  COMMAND $<TARGET_FILE:test_matrix>
  )

add_test(
  NAME test_spglib_utils
  COMMAND $<TARGET_FILE:test_spglib_utils>
  )

add_test(
  NAME test_abf
  COMMAND $<TARGET_FILE:test_abf>
  )

add_test(
  NAME test_rotate
  COMMAND $<TARGET_FILE:test_rotate>
  )

add_test(
  NAME test_mathtools
  COMMAND $<TARGET_FILE:test_mathtools>
  )

add_test(
  NAME test_kmesh
  COMMAND $<TARGET_FILE:test_kmesh>
  )

add_test(
  NAME test_vec
  COMMAND $<TARGET_FILE:test_vec>
  )

add_test(
  NAME test_cell
  COMMAND $<TARGET_FILE:test_cell>
  )

add_test(
  NAME test_trans
  COMMAND $<TARGET_FILE:test_trans>
  )

add_test(
  NAME test_io
  COMMAND $<TARGET_FILE:test_io>
  )

