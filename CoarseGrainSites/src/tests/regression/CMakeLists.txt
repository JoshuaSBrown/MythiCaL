
foreach(PROG test_ToF
    test_1Dsystem)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(regression_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(regression_${PROG} mythical)
endforeach(PROG)

configure_file(test_regression_script.sh test_regression_script.sh COPYONLY)

add_test(test_regression bash test_regression_script.sh)
set_tests_properties(test_regression PROPERTIES LABELS "mythical")
