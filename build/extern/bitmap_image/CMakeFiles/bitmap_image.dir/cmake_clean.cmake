file(REMOVE_RECURSE
  "libbitmap_image.a"
  "libbitmap_image.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang )
  include(CMakeFiles/bitmap_image.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
