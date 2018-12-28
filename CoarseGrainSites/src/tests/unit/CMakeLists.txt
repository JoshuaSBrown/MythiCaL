
foreach(PROG 
    test_identity 
    test_kmc_basin_explorer
    test_kmc_cluster 
    test_kmc_cluster_container
    test_kmc_coarsegrainsystem
    test_kmc_coarsegrainsystem2
    test_kmc_graph_library_adapter
    test_kmc_walker
    test_kmc_rate_container
    test_kmc_site
    test_kmc_site_container)

  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} kmccoarsegrain)
  add_test(unit_${PROG} unit_${PROG})
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "kmccoarsegrain")
endforeach(PROG)
