file(REMOVE_RECURSE
  "libkmc_course_grain.pdb"
  "libkmc_course_grain.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/kmc_course_grain.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
