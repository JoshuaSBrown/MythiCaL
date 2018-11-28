
foreach(PROG test_ToF)

  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(regression_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(regression_${PROG} kmc_course_grain ugly)
  add_test(regression_${PROG} regression_${PROG})
  set_tests_properties(regression_${PROG} PROPERTIES LABELS "kmccoursegrain")
endforeach(PROG)
