#include "ultimaille/attributes.h"
#include "ultimaille/io/by_extension.h"
#include <list>
#include <ultimaille/all.h>

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
        Halfedge nextHalfedge = currentHalfedge.opposite().next();
      
        if (nextHalfedge.active()){
            currentHalfedge = nextHalfedge;

            if (currentHalfedge == startHalfedge)
                return valence;
        }
        else {
            return valence;
        }
    }

    return 4;
}

inline bool isNewDefect(Vertex v, std::vector<int>& defects) {
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

inline int checkTopologicalDisk(FacetAttribute<int>& fa, Quads& m, std::list<int>& patch){
    for (int i : patch)
        fa[Halfedge(m, i).facet()] = 3;

    // Veryfing that we have a topological disk, i.e. all the facets surrounding the inside of the patch are in the patch
    for (int i = 0; i < m.nfacets(); i++){
        if (fa[i] > 0 && fa[i] < 3){ // facet in the patch but not on the outline
            for (Halfedge he : Facet(m, i).iter_halfedges()){
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
    return 1;
}

inline int postPatch(FacetAttribute<int>& fa, Quads& m, std::list<int>& patch, std::list<int>& patchConvexity, CornerAttribute<int>& ca){
    if (checkTopologicalDisk(fa, m, patch) == -1)
        return -1;

    // check if no hard edge has been violated    
    for (Halfedge he : m.iter_halfedges()){
        if (ca[he] == 1 && fa[he.facet()] > 1 && fa[he.opposite().facet()] > 1){
            return -1;
        }
    }

    patchRotationRightToEdge(patch, patchConvexity);

    int nbEdge = 0;
    for (int convexity : patchConvexity){
        if (convexity >= 1)
            nbEdge++;
    }
    
    return nbEdge;
}

inline int bfs(int startFacet, FacetAttribute<int>& fa, Quads& m, CornerAttribute<int>& ca){
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
        if (ca[facetHalfedge] == 1)
            continue;

        // explore each vertex of the facet        
        for (int i = 0; i < 4; i++){
            facetHalfedge = facetHalfedge.next();

            exploringHalfedge = facetHalfedge;

            if (isNewDefect(exploringHalfedge.from(), defectVertices)){
                defectCount++;

                // marking the facets surrounding the last defect vertex
                if (defectCount == 3){
                    for (int i = 0; i < MAX_VALENCE; i++){
                        fa[exploringHalfedge.facet()] = 1;
                        if (ca[exploringHalfedge] == 1)
                            break;
                        exploringHalfedge = exploringHalfedge.opposite().next();
                        
                    }
                    break;
                }
            }
            
            // explore the facets surrounding a given vertex to mark them
            for (int i = 0; i < MAX_VALENCE; i++){ 
                fa[exploringHalfedge.facet()] = 1;
                otherHalfedge = exploringHalfedge.opposite().next();
                if (ca[exploringHalfedge] == 1)
                    break;
                
                if (fa[otherHalfedge.facet()] == 0){
                    facetQueue.push(otherHalfedge.facet());
                    exploringHalfedge = otherHalfedge;
                }
            }
        }
    }

    // find a halfedge in the border of the patch, inside the patch
    int max_iter = 500;
    while(fa[facetHalfedge.facet()] >= 1 && max_iter > 0){
        facetHalfedge = facetHalfedge.next().next().opposite();
        max_iter--;
    }

    facetHalfedge = facetHalfedge.opposite();

    assert(fa[facetHalfedge.facet()]>0);
    assert(fa[facetHalfedge.opposite().facet()]<1);

    return facetHalfedge;
}


inline int getPatch(Halfedge boundaryHe, FacetAttribute<int>& fa, std::list<int>& halfedgePatch, std::list<int>& patchConvexity){
    // We want a list of all the halfedge on the boundary of the patch (information is in fa)
    // boundaryHe is a halfedge on the patch. We start from here and do a rotation outward of the patch to find the next halfedge, and so on until coming back to the start

    assert(boundaryHe >= 0);
    halfedgePatch.clear();
    patchConvexity.clear();

    halfedgePatch.push_back(boundaryHe);
    Halfedge startHalfedge = boundaryHe;
    bool firstIteration = true;


    while (boundaryHe != startHalfedge || firstIteration){
        assert(boundaryHe.opposite().prev().opposite() != -1);
 
        boundaryHe = boundaryHe.opposite();
        for (int i=0; i<MAX_VALENCE; i++){
            boundaryHe = boundaryHe.prev().opposite();
            if (fa[boundaryHe.facet()] >= 1){
                halfedgePatch.push_back(boundaryHe);
                patchConvexity.push_back(i-1);
                break;
            }
        }

        firstIteration = false;
    }  

    halfedgePatch.pop_front();
    return 1;
}
inline int updateBoundaryHe(int& boundaryHe, Halfedge& he, Quads& m, FacetAttribute<int>& fa){
    // Makes sure we have a halfedge on the boundary of the patch

    he = he.next();
    if (fa[he.opposite().facet()] < 1){
        boundaryHe = he;
    } else if (fa[he.next().opposite().facet()] < 1){
        boundaryHe = he.next();
    } else if (fa[he.next().next().opposite().facet()] < 1){
        boundaryHe = he.next().next();
    } else {
        return -1;
    }

    assert(fa[Halfedge(m, boundaryHe).facet()]>0);
    assert(fa[Halfedge(m, boundaryHe).opposite().facet()]<1);
    return 1;
}

inline int makePatchConcave(int& boundaryHe, std::list<int>& patch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, Quads& m, CornerAttribute<int>& ca){
    int max_iter = 100;
    bool hasConcave = true;
    while (hasConcave && max_iter > 0){
        hasConcave = false;
        Halfedge he = Halfedge(m, boundaryHe);
        for (auto [a, b] : zip(patch, patchConvexity)) {
            if (b < 0){
                if (ca[he] == 1)
                    continue;
                he = Halfedge(m, a).opposite();

                fa[he.facet()] = 2;
                hasConcave = true;
            }
        }

        if (updateBoundaryHe(boundaryHe, he, m, fa) == -1)
            return -1;

        getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity);
    }
    return 1;
}

inline int initialPatchConstruction(Vertex v, FacetAttribute<int>& fa, std::list<int>& patch, std::list<int>& patchConvexity, Quads& m, CornerAttribute<int>& ca){
    // constructing a patch with 3 defects with breath-first search
    int boundaryHe = bfs(v.halfedge().facet(), fa, m, ca);
    getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity);
    makePatchConcave(boundaryHe, patch, patchConvexity, fa, m, ca);
    return postPatch(fa, m, patch, patchConvexity, ca);
}

inline int expandPatch(std::list<int>& patch, FacetAttribute<int>& fa, Quads& m, std::list<int>& patchConvexity, CornerAttribute<int>& ca){
    Halfedge he = Halfedge(m, 1);
    for (int i : patch) {
        if (ca[Halfedge(m, i)] == 1)
            continue;
        he = Halfedge(m, i).opposite();
        fa[he.facet()] = 2;
    }

    int boundaryHe = 0;
    if (updateBoundaryHe(boundaryHe, he, m, fa)==-1)
        return -1;

    getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity);

    if (makePatchConcave(he, patch, patchConvexity, fa, m, ca)==-1)
        return -1;
    

    return postPatch(fa, m, patch, patchConvexity, ca);
}
