# Implementing Bunin's algorithm with Ultimaille

Improve quad mesh topology with non-local topological clean-up.

## Compilation

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cd build
make -j
```

- add `#include <cstdint>` on ultimaille/io/geogram.cpp if compilation error is found
- `cmake -B build -DCMAKE_BUILD_TYPE=Debug` for debugging (delete the build folder before)

## Usage

```
./src/main {input path}
graphite output.geogram ../geogram.lua
```

## Roadmap 

- Expanding the patch if remeshing fail (recommended max 5000 facets, or less if working on curvy 3d meshes, cf. Jaal)
- Working with 3 and 4 sided loops
- Using Djisktra to find the second defect
- Explore different seed chosing methods to start the algorithm, and between each iteration
- Smoothing (Laplacian, Quasi-Newton non-linear, other?)
- 3d geometry considerations

## Notes

Profiling :
```
valgrind --tool=callgrind ./src/main
kcachegrind callgrind.out.*
```

## References
DOI:10.1007/978-3-540-34958-7_1
DOI:10.1007/978-3-642-24734-7-28  with source code https://ftp.mcs.anl.gov/pub/fathom/meshkit-docs//index.html
DOI:10.1016/j.proeng.2015.10.137