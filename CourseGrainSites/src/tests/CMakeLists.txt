
foreach(PROG test_kmc_cluster test_kmc_coursegrainsystem test_identity test_kmc_particle test_kmc_site )
  file(GLOB ${PROG}_SOURCES ${PROG}*.cpp)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} kmc_course_grain ugly)
  add_test(unit_${PROG} unit_${PROG})
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "kmccoursegrain")
endforeach(PROG)
