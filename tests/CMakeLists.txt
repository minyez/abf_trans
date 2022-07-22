add_executable(test_sph test_sph.cpp)
add_executable(test_linalg test_linalg.cpp)
add_executable(test_matrix test_matrix.cpp)
add_executable(test_spglib_utils ./test_spglib_utils.cpp)

target_sources(test_sph PRIVATE ../src/sph.cpp)
target_sources(test_spglib_utils PRIVATE ../src/spglib_utils.cpp)

target_link_libraries(test_linalg "${LAPACK_LIBRARIES}")

target_link_libraries(test_matrix "${LAPACK_LIBRARIES}")

target_link_libraries(test_spglib_utils "${SPGLIB_LIBRARIES} ${LAPACK_LIBRARIES}")

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
