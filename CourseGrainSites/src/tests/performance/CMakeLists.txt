
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

foreach(PROG test_kmc_coursegrainsystem)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(performance_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(performance_${PROG} kmc_course_grain ugly)
  add_test(performance_${PROG} performance_${PROG})
  set_tests_properties(performance_${PROG} PROPERTIES LABELS "kmccoursegrain")
endforeach(PROG)

foreach(PROG  
    test_crude_vs_coursegrain
    test_crude_vs_coursegrain_correlated)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(performance_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(performance_${PROG} kmc_course_grain ugly)
  add_test(performance_${PROG} performance_${PROG} 0.08 10 10 10 3 1000 0.001)
  set_tests_properties(performance_${PROG} PROPERTIES LABELS "kmccoursegrain")
endforeach(PROG)

configure_file(test_script_crude_vs_course.sh test_script_crude_vs_course.sh COPYONLY)

