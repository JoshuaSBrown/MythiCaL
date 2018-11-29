
foreach(PROG test_ToF)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(regression_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(regression_${PROG} kmc_course_grain ugly)
endforeach(PROG)

add_test(test_regression test_regression_script.sh)
set_tests_properties(test_regression PROPERTIES LABELS "kmccoursegrain")
configure_file(test_regression_script.sh test_regression_script.sh COPYONLY)
