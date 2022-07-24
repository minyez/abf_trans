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
target_sources(test_rotate PRIVATE ../src/rotate.cpp ../src/constants.cpp)
target_link_libraries(test_rotate "${LAPACK_LIBRARIES}")

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
