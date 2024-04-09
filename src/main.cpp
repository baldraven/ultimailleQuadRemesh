#include "ultimaille/attributes.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <ultimaille/all.h>
#include <list>
#include "patchFinding.h"
#include "remeshing.h"

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;


void animate(Quads& m, int i){
    std::string number = std::to_string(i);
    if (i < 10){
        number = "0" + number;
    }
    std::string s = "../animation/output" + number + ".geogram";
    write_by_extension(s, m);
}

int countDefect(Quads& m){
    int count = 0;
    for (Vertex v: m.iter_vertices()){
        if (getValence(v) != 4){
            count++;
        }
    }
    return count;
}


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

    int defectCountBefore = countDefect(m);

    FacetAttribute<int> fa(m, 0);
    
    // iterating through the mesh until no remesh can be done
    for (int i=0; i < 200; i++){
        bool hasRemeshed = false;

        // iterating through the vertices until finding a defect 
        for (Vertex v: m.iter_vertices()){
            // resetting the attribute
            fa.fill(0);

            if (v.on_boundary())
                continue;
            if (getValence(v) == 4)
                continue;
            
            // constructing the patch and coloring the facets inside
            int boundaryHe = bfs(v.halfedge().facet(), fa, m);
            if (boundaryHe == -1)
                continue;

            std::list<int> patch; // we could resort to Dijsktra to make the patch smaller
            std::list<int> patchConvexity;

            int facetCount = 0;
            int iter = 0;

            // Magic numbers, to tweak
            while (facetCount < 500 && iter < 50){ 

                int edgeCount = completingPatch(boundaryHe, fa, m, patch, patchConvexity);

                facetCount = countFacetsInsidePatch(fa, m.nfacets()); // we could increment the count instead of counting the facets each time

                // remeshing the patch
                 if (edgeCount == 3){

                    if(remeshing3patch(patch, patchConvexity, m, fa, v)){
                        hasRemeshed = true;
                        break;
                    }
                } 

                if (edgeCount == 5){

                    if(remeshing5patch(patch, patchConvexity, m, fa, v)){
                        hasRemeshed = true;
                        break;
                    }
                }



                iter++;
            }
            if (hasRemeshed){
                break;
            }
        }

        animate(m, i);

        if (!hasRemeshed){
            std::cout << "No more defects found after " << i+1 << " remeshing." << std::endl;
            break;
        }

        //debug at the first remesh
        //break;
    }

    write_by_extension("output.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});

    int defectCountAfter = countDefect(m);
    std::cout << "Number of corrected defects: " << defectCountBefore-defectCountAfter << " out of " << defectCountBefore << std::endl;

    return 0;
}


