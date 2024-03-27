#include "helpers.h"
#include "ultimaille/attributes.h"
#include <iostream>
#include <ultimaille/all.h>
#include <list>
#include "patchFinding.h"
#include "remeshing.h"

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

int main() {
    std::string path = getAssetPath();
    Quads m;
    //read_by_extension(path + "catorus_voxel_bound_smooth.geogram", m);
    read_by_extension(path + "cowhead.geogram", m);
    m.connect();

    FacetAttribute<int> fa(m, 0);

    // iterating through the vertices until finding a defect 
    for (int _=0; _ < 1; _++){
        std::cout << "New iteration" << std::endl;

        for (Vertex v: m.iter_vertices()){
            // reset the attribute
            fa.fill(0);

            if (v.on_boundary())
                continue;

            if (getValence(v) == 4)
                continue;
            
            // constructing the patch and coloring the facets inside
            int boundaryHe = bfs(v.halfedge().facet(), fa, m);
            if (boundaryHe == -1)
                continue;


/*             if (v==1025)
                break;  */
 
 
            // expanding the patch to include concave facets. patch is a list of halfedges in the boundary of the patch
            std::list<int> patch;
            std::list<int> patchConvexity;
            int edgeCount = completingPatch(boundaryHe, fa, m, patch, patchConvexity);
            if (edgeCount == 1000)
                break;

            // remeshing the patch
            if (edgeCount == 5){
                remeshing5patch(patch, patchConvexity, m, fa, v);
            }

        }
    }

    write_by_extension("bunnin.geogram", m, {{}, {{"patch", fa.ptr}}, {}});
    return 0;
}
    