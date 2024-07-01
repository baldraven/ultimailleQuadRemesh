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
./{executable path} model={model path} {paramType=paramValue}
```

Optional parameters are :
- *string* **result_path** : sets output meshes path (defauls to *output/*)
- *bool* **animate** : sets whether to export the mesh after each iteration, in a output/animation folder (defaults to *false*)
- *int* **maxPatchSize** : sets the maximum number of facets in a patch to remesh. Higher usually eliminate more defects, but can be slower (defaults to *500*)
- *bool* **cad_mode** : enable a mode that preserve the edges of the mesh (default to *false*)

Alternatively, it can be run from Graphite with [graphite addon loader](https://github.com/ultimaille/graphite-addon-loader).

## How does it work

Starting by the first point of the mesh, it looks for patches with at least 3 singularities (points that have a number of incident edges different from 4), that have either 3, 4 or 5 sides, and remesh them with a single singularity, thanks to [these equations](src/matrixEquations.h).

If it fails, it expands the patch until reaching a maximum size. Then it continues iterating on the other points. If a remesh is done, it restarts iterating from the first point again, continuing until every point is covered and no remesh has been made.

## References
- DOI:10.1007/978-3-540-34958-7_1
- DOI:10.1007/978-3-642-24734-7-28  with source code https://ftp.mcs.anl.gov/pub/fathom/meshkit-docs//index.html
- DOI:10.1016/j.proeng.2015.10.137