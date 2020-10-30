
if(CODE_COVERAGE)
  if( NOT COVERAGE_NAME)
    SET(COVERAGE_NAME "coverage_reports")
  endif()
  if( NOT COVERAGE_PATH)
    SET(COVERAGE_PATH "${CMAKE_BINARY_DIR}/coverage")
  endif()
  SET(GCC_COVERAGE_COMPILE_FLAGS "--coverage")
  SET(GCC_COVERAGE_LINK_FLAGS    "--coverage")
endif()

function(add_coverage_label tests )
  if( CODE_COVERAGE )
    foreach( CHECK_TEST ${tests})
      set(exclude FALSE)
      get_test_property(${CHECK_TEST} LABELS TEST_LABELS)
      foreach( exclude_if_contains_this_label ${ARGN})
        if( ${exclude_if_contains_this_label} IN_LIST TEST_LABELS)
          set(exclude TRUE) 
          continue()
        endif()
      endforeach()
      if(${exclude})
        continue()
      endif()
      set_property(TEST "${CHECK_TEST}" APPEND PROPERTY LABELS "coverage")
    endforeach()
  endif()
endfunction()

function(create_coverage_targets)
  if(CODE_COVERAGE)

    find_program( PATH_GCOV gcov )
    find_program( PATH_LCOV lcov )

    if(NOT PATH_GCOV)
      message(FATAL_ERROR "Unable to build with code coverage gcov was not found.")
    endif() 

    if(NOT PATH_LCOV)
      message(FATAL_ERROR "Unable to build with code coverage lcov was not found.")
    endif() 

    if(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
      message(WARNING "Using code coverage with an optimized build is discouraged, as it may lead to misleading results.")
    endif() 

    message(STATUS "Coverage reports will be placed in ${COVERAGE_PATH}/${COVERAGE_NAME}")

    get_target_property(MYTHICAL_SOURCES mythical SOURCES)
    get_target_property(UNIT_TEST_SOURCES unit_tests SOURCES)

    add_custom_target(coverage)
    add_custom_command(TARGET coverage

      COMMAND echo "====================== Code Coverage ======================"
      COMMAND mkdir -p ${COVERAGE_PATH}/${COVERAGE_NAME}
      COMMAND ${PATH_LCOV} --version
      # Clean
      COMMAND ${PATH_LCOV} --gcov-tool ${PATH_GCOV} --directory ${PROJECT_BINARY_DIR} -b ${PROJECT_SOURCE_DIR} --zerocounters
      # Base report
      COMMAND ctest -L coverage --verbose
      COMMAND ${PATH_LCOV} --gcov-tool ${PATH_GCOV} -c -i -d ${PROJECT_BINARY_DIR} -b ${PROJECT_SOURCE_DIR} -o ${COVERAGE_PATH}/${COVERAGE_NAME}/report.base.old
      # Remove Kokkos and tst info from code coverage
      COMMAND ${PATH_LCOV} --remove ${COVERAGE_PATH}/${COVERAGE_NAME}/report.base.old "*/UGLY/*" "*/tst/*" -o ${COVERAGE_PATH}/${COVERAGE_NAME}/report.base

      # Capture information from test runs
      COMMAND ${PATH_LCOV} --gcov-tool ${PATH_GCOV} --directory ${PROJECT_BINARY_DIR} -b ${PROJECT_SOURCE_DIR} --capture --output-file ${COVERAGE_PATH}/${COVERAGE_NAME}/report.test.old
      # Remove Kokkos and tst info from code coverage
      COMMAND ${PATH_LCOV} --remove ${COVERAGE_PATH}/${COVERAGE_NAME}/report.test.old "*/UGLY/*" "*/tst/*" -o ${COVERAGE_PATH}/${COVERAGE_NAME}/report.test

      # Combining base line counters with counters from running tests
      COMMAND ${PATH_LCOV} --gcov-tool ${PATH_GCOV} -a ${COVERAGE_PATH}/${COVERAGE_NAME}/report.base -a ${COVERAGE_PATH}/${COVERAGE_NAME}/report.test --output-file ${COVERAGE_PATH}/${COVERAGE_NAME}/report.all
      # Remove unneeded reports
      COMMAND rm ${COVERAGE_PATH}/${COVERAGE_NAME}/report.test.old
      COMMAND rm ${COVERAGE_PATH}/${COVERAGE_NAME}/report.base.old
      COMMAND rm ${COVERAGE_PATH}/${COVERAGE_NAME}/report.base
      COMMAND rm ${COVERAGE_PATH}/${COVERAGE_NAME}/report.test
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
      )

    set(UPLOAD_COMMAND "bash <\(curl -s https://codecov.io/bash\) \|\| echo \"code coverage failed to upload\"")
    add_custom_target(coverage-upload) 
    add_custom_command(TARGET coverage-upload
      COMMAND echo "================ Uploading Code Coverage =================="
      # Upload coverage report
      COMMAND ${PROJECT_SOURCE_DIR}/scripts/combine_coverage.sh ${PATH_LCOV} ${PATH_GCOV} ${COVERAGE_PATH}
      COMMAND curl -s https://codecov.io/bash > ${COVERAGE_PATH}/CombinedCoverage/script.coverage
      COMMAND cat ${COVERAGE_PATH}/CombinedCoverage/script.coverage
      COMMAND cd ${COVERAGE_PATH}/CombinedCoverage && bash ${COVERAGE_PATH}/CombinedCoverage/script.coverage -p ${PROJECT_BINARY_DIR} -s ${COVERAGE_PATH}/CombinedCoverage
      WORKING_DIRECTORY ${COVERAGE_PATH}
      )

    if(ENABLE_UNIT_TESTS)
      add_dependencies(coverage unit_tests)
    endif()

  else()
    add_custom_target(coverage)
    add_custom_command(TARGET coverage
      COMMAND echo "====================== Code Coverage ======================"
      COMMENT "Code coverage has not been enabled"
      )
  endif()
endfunction()
