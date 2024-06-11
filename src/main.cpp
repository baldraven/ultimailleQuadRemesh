#include "ultimaille/attributes.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <ultimaille/all.h>
#include <list>
#include "patchFinding.h"
#include "remeshing.h"
#include <filesystem>
#include <unistd.h>
#include "param_parser.h"
#include "ultimaille/io/by_extension.h"


using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

void animate(Quads& m, int i){
    std::cout << "animate: " << i << " | ";
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

bool loadingInput(Quads& m, std::string path){
    read_by_extension(path, m);

    if (m.nverts() == 0) {
        std::cerr << "Error reading file" << std::endl;
        return false;
    }

    m.connect();
    return true;
}

void mainLoop(Quads& m, BVH& bvh, FacetAttribute<int>& fa){
    int MAX_REMESH = 999;
    for (int i=0; i < MAX_REMESH; i++){
        bool hasRemeshed = false;
        
        // iterating through the vertices until finding a defect 
        for (Vertex v: m.iter_vertices()){
            fa.fill(0);

            if (getValence(v) == 4)
                continue;
            
            std::list<int> patch; 
            std::list<int> patchConvexity;
            int edgeCount = initialPatchConstruction(v, fa, patch, patchConvexity, m);
            if (edgeCount == -1)
                continue;

            // trying to remesh and expanding the patch in case of failure until we reach the maximum patch size
            int facetCount = 0;
            int max_iter = 20;
            int MAX_PATCH_FACET_COUNT = 500;
            while (facetCount < MAX_PATCH_FACET_COUNT && max_iter > 0){
              //  std::cout << "JHUJ" << std::endl;
                
                if ( edgeCount == 4 || edgeCount == 3 || edgeCount == 5){
                    write_by_extension("../output/debugB.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});
                    if(remeshingPatch(patch, patchConvexity, edgeCount, m, fa, v, bvh)){
                        hasRemeshed = true;
                        break;
                    }
                }

                edgeCount = expandPatch(patch, fa, m, patchConvexity);
                if (edgeCount == -1){
           //         std::cout << "-1" << std::endl;
                    
                    break; 
                }

               //animate(m, i);
                //std::cout << "REACHED" << std::endl;
                

                facetCount = countFacetsInsidePatch(fa, m.nfacets());
                max_iter--;
            }

            if (hasRemeshed){
                animate(m, i);
                break;
            }
        }

        if (!hasRemeshed){
            std::cout << "No more defects found after " << i << " remeshing." << std::endl;
            break;
        }
    }
}

int main(int argc, char* argv[]) {
    Parameters params;
    params.help = "This addon delete a facet \n yop";
    // Add program parametersw
    params.add("input", "model", "").description("Model to process");
    params.add("string", "result_path", "").type_of_param("system");
    /* Parse program arguments */
    params.init_from_args(argc, argv);
    // Get parameters
    std::string filename = params["model"];
    std::filesystem::path result_path((std::string)params["result_path"]);

    Quads m;
    if (!loadingInput(m, filename))
        return EXIT_SUCCESS;

    // Constructing structure for projecting the new patches on the original mesh
    Triangles mTri = quand2tri(m);
    BVH bvh(mTri);  

    int defectCountBefore = countDefect(m);

    FacetAttribute<int> fa(m, 0);
    mainLoop(m, bvh, fa);

    // Output model to output directory at working dir if no result_path given
    if (result_path.empty() && !std::filesystem::is_directory("output")) {
        std::filesystem::create_directories("output");
        result_path = "output";
    }

    // Get file name and output path
    std::string file = std::filesystem::path(filename).filename().string();
    std::string out_filename = (result_path / file).string();

    write_by_extension(out_filename, m, {{}, {{"patch", fa.ptr}, }, {}});
    std::cout << "Result exported in " << out_filename << std::endl;
    

    int defectCountAfter = countDefect(m);
    std::cout << "Number of corrected defects: " << defectCountBefore-defectCountAfter << " out of " << defectCountBefore << std::endl;

    return 0;
}
