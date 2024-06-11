![build, test](https://github.com/baldraven/ultimailleQuadRemesh/actions/workflows/continuous.yml/badge.svg)

<p align="center">
  <img src="https://i.imgur.com/H9Hh1gZ.jpg"/>
</p>

# Implementing Bunin's algorithm with Ultimaille

Improve quad mesh topology with non-local topological clean-up.

## Compilation

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cd build
make -j
```
## Usage

```
./{executable path} model={model path} result_path={result path}
```

## Roadmap 

- Add integration to Graphite
- Add parameters to the executable, like max patch size, first seed selection or enabling animated output
- Implement edge flipping
- Improve result for CAD-like input meshes, by avoiding remeshing hardedges
- Improve geometry
- Improve patch construction algorithms


## How does it work

Starting by the first point of the mesh, it looks for patches with at least 3 singularities (points that have a number of incident edges different from 4), that have either 3, 4 or 5 sides, and try to remesh by finding a new connectivity inside the patch, thanks to [these equations](src/matrixEquations.h).

If it fails, it expands the patch until reaching a maximum size (by default 500 facets). Then it continues iterating on the other points. If a remesh is done, it restarts iterating from the first point again, continuing until every point is covered and no remesh has been made.

## References
- DOI:10.1007/978-3-540-34958-7_1
- DOI:10.1007/978-3-642-24734-7-28  with source code https://ftp.mcs.anl.gov/pub/fathom/meshkit-docs//index.html
- DOI:10.1016/j.proeng.2015.10.137