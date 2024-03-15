file(REMOVE_RECURSE
  "libpcg32.a"
  "libpcg32.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/pcg32.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
