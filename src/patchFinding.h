#include "ultimaille/io/by_extension.h"
#include <list>
#include <ultimaille/all.h>
#include <unistd.h>

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

const int MAX_VALENCE = 8;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// patch finding

inline int getValence(Vertex vertex){
    int valence = 0;
    Halfedge currentHalfedge = vertex.halfedge();
    Halfedge startHalfedge = currentHalfedge;

    for (int i = 0; i < MAX_VALENCE; i++){
        valence++; 
        Halfedge nextHalfedge = currentHalfedge;
      
        if (nextHalfedge.opposite() == -1)
            return 4;

        nextHalfedge = nextHalfedge.opposite().next();
        currentHalfedge = nextHalfedge;

        if (currentHalfedge == startHalfedge)
            return valence;
    }

    assert(false);    
    return 4;
}

inline bool isDefect(Vertex v, std::vector<int>& defects) {
    if (getValence(v) != 4 && 
        std::find(defects.begin(), defects.end(), v) == defects.end()) {
        defects.push_back(v);
        return true;
    }
    return false;
}

inline void patchRotationRightToEdge(std::list<int>& patch, std::list<int>& patchConvexity){
    for (int i=0; i<int(patchConvexity.size()); i++){
        if (patchConvexity.front() < 1){
            patch.splice(patch.begin(), patch, std::prev(patch.end()));
            patchConvexity.splice(patchConvexity.begin(), patchConvexity, std::prev(patchConvexity.end()));
        } else {
            break;
        }
    }
}

inline int countFacetsInsidePatch(FacetAttribute<int>& fa, int nfacets){
    int nbFacetInsidePatch = 0;
    for (int f = 0; f < nfacets; f++){
        if (fa[f] > 1)
            nbFacetInsidePatch++;
    }
    return nbFacetInsidePatch;
}

inline int postPatch(FacetAttribute<int>& fa, Quads& m, std::list<int>& patch, std::list<int>& patchConvexity){
    for (int i : patch)
        fa[Halfedge(m, i).facet()] = 3;

    // Veryfing that we have a topological disk, i.e. all the facets surrounding the inside of the patch are in the patch
    for (int i = 0; i < m.nfacets(); i++){
        if (fa[i] > 0 && fa[i] < 3){ // facet in the patch but not on the outline
            for (Halfedge he : Facet(m, i).iter_halfedges()){
                if (he.opposite() != -1)
                    continue;
                if (fa[he.opposite().facet()] < 1){
                    return -1;
                }
            }
        }
    }

    for (int i : patch){
        if (fa[Halfedge(m, i).facet()] == 3)
            fa[Halfedge(m, i).facet()] = 2;
    }

    patchRotationRightToEdge(patch, patchConvexity);

    write_by_extension("../output/debugS.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});

    int nbEdge = 0;
    for (int convexity : patchConvexity){
        if (convexity >= 1)
            nbEdge++;
    }
    return nbEdge;
}

inline bool onBoundaryOfMeshOrPatch(Halfedge he, FacetAttribute<int>& fa){
    return he.opposite() == -1 || fa[he.opposite().facet()] < 1;
}

inline int bfs(int startFacet, FacetAttribute<int>& fa, Quads& m){
    std::queue<int> facetQueue;
    int defectCount = 0;
    std::vector<int> defectVertices;

    facetQueue.push(startFacet);
    Halfedge facetHalfedge = Facet(m, startFacet).halfedge();
    Halfedge exploringHalfedge = facetHalfedge;
    Halfedge otherHalfedge = exploringHalfedge;

    while (!facetQueue.empty() && defectCount < 3){
        facetHalfedge = Facet(m, facetQueue.front()).halfedge().prev();
        facetQueue.pop();

        // explore each vertex of the facet        
        for (int i = 0; i < 4; i++){
            facetHalfedge = facetHalfedge.next();
            exploringHalfedge = facetHalfedge;

            if (isDefect(exploringHalfedge.from(), defectVertices)){
                defectCount++;

                // marking the facets surrounding the last defect vertex
                if (defectCount == 3){
                    for (int i = 0; i < MAX_VALENCE; i++){
                        fa[exploringHalfedge.facet()] = 1;
                        if (exploringHalfedge.opposite() != -1)
                            exploringHalfedge = exploringHalfedge.opposite().next();
                    }
                    break;
                }
            }
            
            // explore the facets surrounding a given vertex to mark them
            for (int i = 0; i < MAX_VALENCE; i++){ 
                fa[exploringHalfedge.facet()] = 1;
                
                if (exploringHalfedge.opposite() == -1)
                    break;

                otherHalfedge = exploringHalfedge.opposite().next();
                
                if (otherHalfedge.active() && fa[otherHalfedge.facet()] == 0){
                    facetQueue.push(otherHalfedge.facet());
                exploringHalfedge = otherHalfedge;
                }
            }
        }
    }

    // find a halfedge in the border of the patch, inside the patch
    int max_iter = 500;
    while(!onBoundaryOfMeshOrPatch(facetHalfedge, fa) && max_iter > 0){
        facetHalfedge = facetHalfedge.opposite().next().next();
        max_iter--;
    }

    assert(max_iter > 0);
    assert(facetHalfedge != -1);

    //write_by_extension("../outputL.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});
    //std::cout << facetHalfedge.from() << " | " << facetHalfedge.to() << std::endl;

    
    assert(fa[facetHalfedge.facet()]>0);
    assert(onBoundaryOfMeshOrPatch(facetHalfedge, fa));

    return facetHalfedge;
}

inline int getPatch(Halfedge boundaryHe, FacetAttribute<int>& fa, std::list<int>& halfedgePatch, std::list<int>& patchConvexity, Quads& m){
    // We want a list of all the halfedge on the boundary of the patch (information is in fa)
    // boundaryHe is a halfedge on the patch. We start from here and do a rotation outward of the patch to find the next halfedge, and so on until coming back to the start

    //std::cout << boundaryHe.from() << boundaryHe.to() << std::endl;
    

    assert(boundaryHe >= 0);
    halfedgePatch.clear();
    patchConvexity.clear();

    halfedgePatch.push_back(boundaryHe);
    Halfedge startHalfedge = boundaryHe;
    bool firstIteration = true;
    int max_iteration = 500;

    //std::cout << "enter" << std::endl;
    //std::cout << boundaryHe.from() <<" | " << boundaryHe.to() << std::endl;

    while ((boundaryHe != startHalfedge || firstIteration) && max_iteration > 0){
        firstIteration = false;
        max_iteration--;
        //std::cout << boundaryHe.from() << " " << boundaryHe.to() << std::endl;
    
        boundaryHe = boundaryHe.next();
        if (onBoundaryOfMeshOrPatch(boundaryHe, fa)){
                //std::cout << "Break: " << boundaryHe.from() <<" | " << boundaryHe.to() << std::endl;
                halfedgePatch.push_back(boundaryHe);
                patchConvexity.push_back(-1); // TODO : check this
                continue;
        }
        for (int i=0; i<MAX_VALENCE; i++){
       //std::cout << boundaryHe.from() <<" | " << boundaryHe.to() << std::endl;

            // TODO: factorise
            if (boundaryHe.opposite() == -1){
                halfedgePatch.push_back(boundaryHe);
                patchConvexity.push_back(-1);
                break;
            }

            boundaryHe = boundaryHe.opposite().next();
            if (onBoundaryOfMeshOrPatch(boundaryHe, fa)){
               // std::cout << "enter2" << std::endl; 
                halfedgePatch.push_back(boundaryHe);
                patchConvexity.push_back(i-1);
                break;
            }
        }
    }  
    assert(max_iteration > 0);
    //write_by_extension("../animation/debugZ.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});

/*     std::cout << std::endl;
    
    for (int i : halfedgePatch){
        std::cout << Halfedge(m, i).from() << " " << Halfedge(m, i).to() << std::endl;
    }
     */
    
    for (int i : halfedgePatch){
        assert(fa[Halfedge(m, i).facet()] >= 1);
        assert(onBoundaryOfMeshOrPatch(Halfedge(m, i), fa)); // TODO: assert for borders
    }

    assert(halfedgePatch.front()==halfedgePatch.back());
    halfedgePatch.pop_front();
    return 1;
}


inline int updateBoundaryHeNonLinear(int& boundaryHe, Halfedge& he, Quads& m, FacetAttribute<int>& fa){
    // Makes sure we have a halfedge on the boundary of the patch
    //std::cout << he.from() << " | " << he.to() << std::endl;

    //he = he.next();
    //std::cout << "HA: " << he.facet() <<  std::endl;


    if (he == -1 || fa[he.facet()] < 1){
        boundaryHe = he.opposite();
    }  else if (he.opposite() == -1 || fa[he.opposite().facet()] < 1){
        boundaryHe = he;
    } else if (he.next().opposite() == -1 || fa[he.next().opposite().facet()] < 1){
        boundaryHe = he.next();
    } else if (he.next().next().opposite() == -1 || fa[he.next().next().opposite().facet()] < 1){
        boundaryHe = he.next().next();
    } else if (he.next().next().next().opposite() == -1 || fa[he.next().next().next().opposite().facet()] < 1){
        boundaryHe = he.next().next().next();
    } else {
        std::cout << he.from() << " | " << he.to() << " | " << he.facet() << std::endl;
        write_by_extension("../output/debugF.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});
        return -1;
    }
    //std::cout << Halfedge(m, boundaryHe).from() << " | " << Halfedge(m, boundaryHe).to() << std::endl;
   // write_by_extension("../animation/debugU.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});
    
    assert(fa[Halfedge(m, boundaryHe).facet()]>0);
    assert(onBoundaryOfMeshOrPatch(Halfedge(m, boundaryHe), fa));
    return 1;
}

inline int updateBoundaryHe(int& boundaryHe, Halfedge& he, Quads& m, FacetAttribute<int>& fa){
    for (Facet f : m.iter_facets()){
        if (fa[f] >= 1){
            for (Halfedge h : f.iter_halfedges()){
                if (onBoundaryOfMeshOrPatch(h, fa)){
                    boundaryHe = h;
                    return 1;
                }
            }
        }
    }
    write_by_extension("../output/debugG.geogram",m, {{}, {{"patch", fa.ptr}, }, {}});
    return -1;
}

inline int makePatchConcave(int& boundaryHe, std::list<int>& patch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, Quads& m){
    int max_iter = 100;
    bool hasConcave = true;
    write_by_extension("../animation/debugPreConcave.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});

    int patchSize = 0;
    for (Facet f : m.iter_facets()){
        if (fa[f] >= 1)
            patchSize++;
    }

    // TODO: patch Size as a global
    while (hasConcave && max_iter > 0 && patchSize < 500){

        hasConcave = false;
        Halfedge he = Halfedge(m, boundaryHe);
        //std::cout << "HA: " << he.from() << " | " << he.to() << std::endl;
        
        
        for (auto [a, b] : zip(patch, patchConvexity)) {
            he = Halfedge(m, a);
            if (b < 0 && he.opposite() != -1){           
                he = he.opposite();

                // Remesh more but can make some ugly cases, we have to verify that we can fix this with better geometry, or uncomment the code
                if (fa[he.next().next().opposite().facet()] >= 1 || fa[he.facet()] >= 1)
                    return -1;

                if (fa[he.facet()] < 1){
                    fa[he.facet()] = 2;
                    patchSize++;
                }
                

                hasConcave = true;
            }
        }

        if (updateBoundaryHe(boundaryHe, he, m, fa) == -1)
            return -1;

        //write_by_extension("../animation/debug.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});

        getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity, m);
    }
    write_by_extension("../animation/debugPostConcave.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});
    assert(max_iter > 0);
    return 1;
}

inline int initialPatchConstruction(Vertex v, FacetAttribute<int>& fa, std::list<int>& patch, std::list<int>& patchConvexity, Quads& m){
    // constructing a patch with 3 defects with breath-first search
    int boundaryHe = bfs(v.halfedge().facet(), fa, m);
    getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity, m);
    if (makePatchConcave(boundaryHe, patch, patchConvexity, fa, m) == -1)
        return -1;
    return postPatch(fa, m, patch, patchConvexity);
}

inline int expandPatch(std::list<int>& patch, FacetAttribute<int>& fa, Quads& m, std::list<int>& patchConvexity){

    Halfedge he = Halfedge(m, 1);
    for (int i : patch) {
        he = Halfedge(m, i);
        if (he.opposite() != -1)
            fa[he.facet()] = 2;
    }

    int boundaryHe = 0;
    if (updateBoundaryHe(boundaryHe, he, m, fa)==-1)
        return -1;


    getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity, m);

    if (makePatchConcave(he, patch, patchConvexity, fa, m)==-1)
        return -1;
    

    return postPatch(fa, m, patch, patchConvexity);
}
