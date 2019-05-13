
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

configure_file(test_script_crude_vs_coarse.sh test_script_crude_vs_coarse.sh COPYONLY)
configure_file(test_script_compare_cluster_vs_nocluster.sh test_script_compare_cluster_vs_nocluster.sh COPYONLY)

foreach(PROG test_kmc_coarsegrainsystem)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(performance_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(performance_${PROG} kmccoarsegrain)
  add_test(performance_${PROG} performance_${PROG})
  set_tests_properties(performance_${PROG} PROPERTIES LABELS "kmccoarsegrain")
endforeach(PROG)

foreach(PROG  
    test_crude_vs_coarsegrain
    test_crude_vs_coarsegrain_correlated)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(performance_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(performance_${PROG} kmccoarsegrain)
  add_test(performance_${PROG} performance_${PROG} 0.08 10 10 3 1000 0.001)
  set_tests_properties(performance_${PROG} PROPERTIES LABELS "kmccoarsegrain")
endforeach(PROG)

foreach(PROG test_crude_vs_coarsegrain_tof)
	file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(performance_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(performance_${PROG} kmccoarsegrain)
  add_test(performance_${PROG} performance_${PROG} 0.08 10 10 3 1000 0.000001 0.05)
  set_tests_properties(performance_${PROG} PROPERTIES LABELS "kmccoarsegrain")
endforeach(PROG)

foreach(PROG  
    test_coarsegrain_cluster_vs_nocluster_toy
    test_coarsegrain_cluster_vs_nocluster)
  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(performance_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(performance_${PROG} kmccoarsegrain)
  add_test(performance_${PROG} bash test_script_compare_cluster_vs_nocluster.sh)
  set_tests_properties(performance_${PROG} PROPERTIES LABELS "kmccoarsegrain")
endforeach(PROG)



