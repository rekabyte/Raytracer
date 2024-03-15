file(REMOVE_RECURSE
  "liblinalg.a"
  "liblinalg.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/linalg.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
