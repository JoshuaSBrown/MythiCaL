
list( APPEND TEST_SOURCES
    catch_main.cpp
    test_identity.cpp 
    test_basin_explorer.cpp
    test_cluster.cpp 
    test_cluster_container.cpp
    test_coarsegrainsystem.cpp
    test_coarsegrainsystem2.cpp
    test_cubic_lattice.cpp
    test_graph_library_adapter.cpp
    test_queue.cpp
    test_walker.cpp
    test_rate_container.cpp
    test_site.cpp
    test_site_container.cpp)

  add_executable(unit_tests ${TEST_SOURCES})
  target_link_libraries(unit_tests PUBLIC mythical Catch2::Catch2)
  ParseAndAddCatchTests(unit_tests)
