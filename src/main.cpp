#include "ultimaille/attributes.h"
#include <cstdlib>
#include <iostream>
#include <string>
#include <ultimaille/all.h>
#include <list>
#include "patchFinding.h"
#include "remeshing.h"
#include <filesystem>


using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

void animate(Quads& m, int i){
    std::string number = std::to_string(i);
    if (i < 10){
        number = "000" + number;
    }
    else if (i < 100){
        number = "00" + number;
    }
    else if (i < 1000){
        number = "0" + number;
    }
    std::string s = "../animation/output" + number + ".geogram";
    write_by_extension(s, m);
}

inline void animateDebug(Quads& m, int i, FacetAttribute<int>& fa){
    std::string number = std::to_string(i);
    if (i < 10){
        number = "000" + number;
    }
    else if (i < 100){
        number = "00" + number;
    }
    else if (i < 1000){
        number = "0" + number;
    }
    std::string s = "../animation/outputDEBUG" + number + ".geogram";
    write_by_extension(s, m, {{}, {{"patch", fa.ptr}, }, {}});
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

Triangles quand2tri(Quads& m){
    Triangles m2;
    m2.points.create_points(m.nverts());
    for(Vertex v : m.iter_vertices()){
        m2.points[v]= v.pos();
    }
    m2.create_facets(m.nfacets()*2);
    for (auto f: m.iter_facets()){
        m2.vert(2*f, 0) = m.vert(f, 0);
        m2.vert(2*f, 1) = m.vert(f, 1);
        m2.vert(2*f, 2) = m.vert(f, 2);

        m2.vert(2*f+1, 0) = m.vert(f, 0);
        m2.vert(2*f+1, 1) = m.vert(f, 2);
        m2.vert(2*f+1, 2) = m.vert(f, 3);
    }
    return m2;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return EXIT_SUCCESS;
    }

    std::string path = argv[1];
    Quads m;

    read_by_extension(path, m);

    if (m.nverts() == 0) {
        std::cerr << "Error reading file" << std::endl;
        return EXIT_SUCCESS;
    }

    m.connect();

    int defectCountBefore = countDefect(m);


    Triangles mTri = quand2tri(m);
    BVH bvh(mTri);


    FacetAttribute<int> fa(m, 0);
    
    // iterating through the mesh until no remesh can be done
    for (int i=0; i < 999; i++){
        bool hasRemeshed = false;
        
        // We'll try something : instead of having a fixed max patch size, we expand it after each failure until we reach the max size
        int maxPatchSize = 500;

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

            while (facetCount < maxPatchSize && iter < 20){ // TODO: tweak those magic numbers

                int edgeCount = completingPatch(boundaryHe, fa, m, patch, patchConvexity, i, v, iter);
                if (edgeCount == -1){
                    break; 
                }

                facetCount = countFacetsInsidePatch(fa, m.nfacets()); // we could increment the count instead of counting the facets each time


                // remeshing the patch
                if ( edgeCount == 4 || edgeCount == 3 || edgeCount == 5){
                    if(remeshingPatch(patch, patchConvexity, edgeCount, m, fa, v, bvh)){
                        hasRemeshed = true;
                        break;
                    }
                }


                iter++;
            }
            if (hasRemeshed){
                std::cout << "animate: " << i << std::endl;
                
                animate(m, i);
                break;
            }
        }

        if (!hasRemeshed && maxPatchSize >= 500){
            std::cout << "No more defects found after " << i+1 << " remeshing." << std::endl;
            break;
        }
        if(!hasRemeshed)
            maxPatchSize += 100;
    }

    std::filesystem::create_directory("output");
    write_by_extension("output/output.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});

    int defectCountAfter = countDefect(m);
    std::cout << "Number of corrected defects: " << defectCountBefore-defectCountAfter << " out of " << defectCountBefore << std::endl;

    return 0;
}
