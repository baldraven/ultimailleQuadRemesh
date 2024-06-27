#include "ultimaille/attributes.h"
#include <cstdlib>
#include <iostream>
#include <numbers>
#include <string>
#include <ultimaille/all.h>
#include <list>
#include "patchFinding.h"
#include "remeshing.h"
#include <filesystem>
#include "param_parser.h"
#include "ultimaille/primitive_geometry.h"
#include "ultimaille/surface.h"

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

void animate(Quads& m, int i, std::string path){
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
    std::string s = path + "/" + number + ".geogram";
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

void edgeFlipping(Quads& m, CornerAttribute<int>& ca){

    bool hasFlipped = true;
    int max_iter = 20;
    while (hasFlipped && max_iter > 0){
        max_iter--;
        hasFlipped = false;
        for (Halfedge he: m.iter_halfedges()){

            if (ca[he] == 1)
                continue;

            Vertex a = he.from();
            Vertex b = he.to();
            int NEa = getValence(a);
            int NEb = getValence(b);
            if (NEa + NEb >= 9){
                Vertex d = he.next().to();
                Vertex f = he.next().next().to();
                Vertex c = he.opposite().next().to();
                Vertex e = he.opposite().next().next().to();

                int NEc = getValence(c);
                int NEe = getValence(e);
                int NEd = getValence(d);
                int NEf = getValence(f);
                if ( (NEa + NEb) - (NEc + NEd) >= (NEa + NEb) - (NEe + NEf) && (NEa + NEb) - (NEc + NEd) >= 3){ // sure about the parenthesis of the && ? 
                    int facet1 = he.facet();
                    int facet2 = he.opposite().facet();

                    for (Halfedge h: Facet(m, facet1).iter_halfedges())
                        if (ca[h] == 1)
                            continue;
                    for (Halfedge h: Facet(m, facet2).iter_halfedges())
                        if (ca[h] == 1)
                            continue;
                    
                    // New quad creation
                    m.conn->create_facet({c, e, b, d});
                    m.conn->create_facet({c, d, f, a});

                    // Cleanup
                    m.conn.get()->active[facet1] = false;
                    m.conn.get()->active[facet2] = false;

                    hasFlipped = true;
                    m.compact(true);
                }
            }
        }
    }
    assert(max_iter > 0);
}

void markHardEdges(Quads& m, CornerAttribute<int>& hardEdges){
    for (Halfedge he: m.iter_halfedges()){

        if (he.opposite() == -1 || hardEdges[he] == 1){
            hardEdges[he] = 1;
            continue;
        }
        
        vec3 n1 = he.facet().geom<Quad3>().normal();
        vec3 n2 = he.opposite().facet().geom<Quad3>().normal();
        n1.normalize();
        n2.normalize();
        double ndot = (n1.x*n2.x + n1.y*n2.y + n1.z*n2.z); 
        double angle = atan2(cross(n1,n2).norm(), ndot);

        if (angle >= std::numbers::pi/4){
            hardEdges[he] = 1;
            hardEdges[he.opposite()] = 1;
        }
    }
}

void mainLoop(Quads& m, BVH& bvh, FacetAttribute<int>& fa, bool ANIMATE, std::string animationPath, int MAXPATCHSIZE, CornerAttribute<int>& ca, bool CAD_MODE){

    if (CAD_MODE)
        markHardEdges(m, ca);

    edgeFlipping(m, ca);

    int i = 0;
    while(true){
        i++;
        bool hasRemeshed = false;
        
        // iterating through the vertices until finding a defect 
        for (Vertex v: m.iter_vertices()){
            fa.fill(0);

            if (getValence(v) == 4)
                continue;

            if (ca[v.halfedge()] == 1 ||
                ca[v.halfedge().next()] == 1 ||
                ca[v.halfedge().next().next()] == 1 ||
                ca[v.halfedge().next().next().next()] == 1)
                continue;
            
            std::list<int> patch; 
            std::list<int> patchConvexity;
            int edgeCount = initialPatchConstruction(v, fa, patch, patchConvexity, m, ca);
            if (edgeCount == -1)
                continue;

            // trying to remesh and expanding the patch in case of failure, until we reach the maximum patch size
            int facetCount = 0;
            int max_iter = 20;
            while (facetCount < MAXPATCHSIZE && max_iter > 0){
                
                if ( edgeCount == 4 || edgeCount == 3 || edgeCount == 5){
                    if(remeshingPatch(patch, patchConvexity, edgeCount, m, fa, v, bvh)){
                        hasRemeshed = true;
                        break;
                    }
                }

                edgeCount = expandPatch(patch, fa, m, patchConvexity, ca);
                if (edgeCount == -1){
                    break; 
                }

                facetCount = countFacetsInsidePatch(fa, m.nfacets());
                max_iter--;
            }

            if (hasRemeshed){
                if (ANIMATE)
                    animate(m, i, animationPath);
                break;
            }
        }

        if (!hasRemeshed){
            std::cout << "No more valid patch found." << std::endl;
            break;
        }
    }
}

int main(int argc, char* argv[]) {
    Parameters params;
    params.help = "This addon correct defects in a quad mesh.";
    params.add("input", "model", "").description("Model to process");
    params.add("string", "result_path", "").type_of_param("system");
    params.add("bool", "animate", "false").description("Export the mesh after each iteration");
    params.add("int", "maxPatchSize", "500").description("Maximum number of facets in a patch to remesh");
    params.add("bool", "cad_mode", "false").description("Respect the sharp angles of the mesh");
    params.init_from_args(argc, argv);

    std::string filename = params["model"];
    std::filesystem::path result_path((std::string)params["result_path"]);
    bool ANIMATE = params["animate"];
    int MAXPATCHSIZE = params["maxPatchSize"];
    bool CAD_MODE = params["cad_mode"];

    Quads m;
    if (!loadingInput(m, filename))
        return EXIT_SUCCESS;

    int defectCountBefore = countDefect(m);

    FacetAttribute<int> fa(m, 0);
    CornerAttribute<int> hardEdges(m, 0);

    if (result_path.empty() && !std::filesystem::is_directory("output")) {
        std::filesystem::create_directories("output");
        result_path = "output";
    }
    std::string animationPath = result_path.string()+"/animation";

    if (ANIMATE && !std::filesystem::is_directory(animationPath)) {
        std::cerr << "Error: " << result_path << " is not a directory" << std::endl;
        std::filesystem::create_directories(animationPath);
    }

    /////////////////////////////////////////////////////////////////////////////////

    // Constructing structure for projecting the new patches on the original mesh
    Triangles mTri = quand2tri(m);
    BVH bvh(mTri);  
    mainLoop(m, bvh, fa, ANIMATE, animationPath, MAXPATCHSIZE, hardEdges, CAD_MODE);

    /////////////////////////////////////////////////////////////////////////////////

    std::string file = std::filesystem::path(filename).filename().string();
    std::string out_filename = (result_path / file).string();
    write_by_extension(out_filename, m, {{}, {{"patch", fa.ptr}, }, {{"hardedges", hardEdges.ptr},}});
    std::cout << "Result exported in " << out_filename << std::endl;

    int defectCountAfter = countDefect(m);
    int percent = 100*(defectCountBefore-defectCountAfter)/defectCountBefore;
    std::cout << "Number of corrected defects: " << defectCountBefore-defectCountAfter << " out of " << defectCountBefore << " (" << percent << ")" << std::endl;

    return percent;
}
