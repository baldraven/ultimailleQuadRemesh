#include "ultimaille/attributes.h"
#include <cstdlib>
#include <iostream>
#include <ultimaille/all.h>
#include <list>
#include "patchFinding.h"
#include "remeshing.h"

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return EXIT_FAILURE;
    }

    std::string path = argv[1];
    Quads m;

    read_by_extension(path, m);

    if (m.nverts() == 0) {
        std::cerr << "Error reading file" << std::endl;
        return EXIT_FAILURE;
    }

    m.connect();

    int count = 0;
    for (Vertex v: m.iter_vertices()){
        if (getValence(v) != 4){
            count++;
        }
    }

    FacetAttribute<int> fa(m, 0);
    
    // iterating through the mesh until no remesh can be done
    for (int i=0; i < 200; i++){
        bool hasRemeshed = false;

        // iterating through the vertices until finding a defect 
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

            // expanding the patch to include concave facets. patch is a list of halfedges in the boundary of the patch
            std::list<int> patch;
            std::list<int> patchConvexity;
            int edgeCount = completingPatch(boundaryHe, fa, m, patch, patchConvexity);

            // remeshing the patch
            if (edgeCount == 5){
                if(remeshing5patch(patch, patchConvexity, m, fa, v)){
                    hasRemeshed = true;
                    break;
                }
            }
        }

        if (!hasRemeshed){
            std::cout << "No more defects found after " << i+1 << " remeshing." << std::endl;
            break;
        }
    }

    write_by_extension("output.geogram", m, {{}, {{"patch", fa.ptr}}, {}});

    for (Vertex v: m.iter_vertices()){
        if (getValence(v) != 4){
            count--;
        }
    }
    std::cout << "Number of corrected defects: " << count << std::endl;

    return 0;
}
    