
foreach(PROG 
    test_identity 
    test_basin_explorer
    test_cluster 
    test_cluster_container
    test_coarsegrainsystem
    test_coarsegrainsystem2
    test_cubic_lattice
    test_graph_library_adapter
    test_queue
    test_walker
    test_rate_container
    test_site
    test_site_container)

  file(GLOB ${PROG}_SOURCES ${PROG}.cpp)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} mythical Catch2::Catch2)
  add_test(unit_${PROG} unit_${PROG})
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "mythical")
endforeach(PROG)
