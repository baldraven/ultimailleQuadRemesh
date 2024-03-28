Implementing Bunin's algorithm on Ultimaille

## Installation
```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cd build
make -j
./src/bunnin
graphite bunnin.geogram ../geogram.lua
```

add `#include <cstdint>` on ultimaille/io/geogram.cpp if compilation error is found
`cmake -B build -DCMAKE_BUILD_TYPE=Debug` for debugging (delete the build folder before)

## References
DOI:10.1007/978-3-540-34958-7_1
DOI:10.1007/978-3-642-24734-7-28  with source code https://ftp.mcs.anl.gov/pub/fathom/meshkit-docs//index.html
DOI:10.1016/j.proeng.2015.10.137