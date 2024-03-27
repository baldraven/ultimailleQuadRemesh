#include <list>
#include <ultimaille/all.h>

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

const int MAX_VALENCE = 8;
const int FACET_VERTICES = 4;

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
            if (currentHalfedge.from() == currentHalfedge.to()){ 
                return 4;
            }

            currentHalfedge = nextHalfedge;

            if(currentHalfedge == startHalfedge){
                return valence;
            }
        }
    }
    throw std::runtime_error("Error in valence calculation for vertex: " + std::to_string(vertex));
    return 0;
}

inline bool isDefect(Vertex v, std::vector<int>& defects) {
    if (getValence(v) != FACET_VERTICES && 
        std::find(defects.begin(), defects.end(), v) == defects.end()) {
        defects.push_back(v);
        return true;
    }
    return false;
}


// TODO: add the last defect into the patch
inline int bfs(int startFacet, FacetAttribute<int>& facetAttributes, Quads& mesh){
    std::queue<int> facetQueue;
    int defectCount = 0;
    std::vector<int> defectVertices;

    facetQueue.push(startFacet);
    Halfedge facetHalfedge = Facet(mesh, startFacet).halfedge();
    Halfedge exploringHalfedge = facetHalfedge;
    Halfedge otherHalfedge = exploringHalfedge;

    while (!facetQueue.empty() && defectCount < 3){
        facetHalfedge = Facet(mesh, facetQueue.front()).halfedge().prev();
        facetQueue.pop();

        // explore each vertex of the facet        
        for (int i = 0; i < FACET_VERTICES; i++){
            facetHalfedge = facetHalfedge.next();

            if (facetHalfedge.from().on_boundary())
                return -1;

            exploringHalfedge = facetHalfedge;

            if (isDefect(exploringHalfedge.from(), defectVertices)) {
                defectCount++;
            }

            // explore the facets surrounding a given vertex to mark them
            for (int i = 0; i < MAX_VALENCE; i++){
                facetAttributes[exploringHalfedge.facet()] = 1;
                
                otherHalfedge = exploringHalfedge.opposite().next();
                if (otherHalfedge.active() && facetAttributes[otherHalfedge.facet()] == 0){
                    facetQueue.push(otherHalfedge.facet());
                    exploringHalfedge = otherHalfedge;
                }
            }
        }
    }

    // find a halfedge in the border of the patch, inside
    // I THINK WE WILL BE ABLE TO MAKE IT A FUNCTION
    while(facetAttributes[facetHalfedge.facet()] >= 1){
        facetHalfedge = facetHalfedge.next().next().opposite();
    }
    facetHalfedge = facetHalfedge.opposite();

    assert(facetAttributes[facetHalfedge.facet()]>0);
    assert(facetAttributes[facetHalfedge.opposite().facet()]<1);

    return facetHalfedge;
}

 inline int getPatch(Halfedge boundaryHe, FacetAttribute<int>& facetAttributes, std::list<int>& halfedgePatch, std::list<int>& patchConvexity){
    assert(boundaryHe > 0);
    halfedgePatch.clear();
    patchConvexity.clear();

    // cycle through the patch to get all the halfedges inside a double linked list
    halfedgePatch.push_back(boundaryHe);
    Halfedge startHalfedge = boundaryHe;
    bool firstIteration = true;

    while (boundaryHe != startHalfedge || firstIteration){
        // boundary issue
        if (boundaryHe.opposite().prev().opposite() == -1){
            return -1;
        }

        boundaryHe = boundaryHe.opposite();
        for (int i=0; i<MAX_VALENCE; i++){
            boundaryHe = boundaryHe.prev().opposite();
            if (facetAttributes[boundaryHe.facet()] >= 1){
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


inline void addConcaveFaces(int& boundaryHe, std::list<int>& patch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, bool& hasConcave, Quads& m){
    hasConcave = false;
    for (auto [a, b] : zip(patch, patchConvexity)) {
        if (b < 0){           
            Halfedge he = Halfedge(m, a).opposite();
            fa[he.facet()] = 2;


            if (!hasConcave){
                // updating the boundaryHe -- we just have to do it one time actually
                he = he.next();
                if (fa[he.opposite().facet()] < 1){
                    boundaryHe = he;
                } else if (fa[he.next().opposite().facet()] < 1){
                    boundaryHe = he.next();
                } else if (fa[he.next().next().opposite().facet()] < 1){
                    boundaryHe = he.next().next();
                } else {
                    std::cout << "?????" << std::endl;
                    continue;
                }
            }

            hasConcave = true;

        }
    }
    assert(fa[Halfedge(m, boundaryHe).facet()]>0);
    assert(fa[Halfedge(m, boundaryHe).opposite().facet()]<1);

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

inline int completingPatch(int boundaryHe, FacetAttribute<int>& fa, Quads& m, std::list<int>& patch, std::list<int>& patchConvexity){
    bool hasConcave = true;
    while(hasConcave){
        if (getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity) == -1){
            return -1;
        }
        addConcaveFaces(boundaryHe, patch, patchConvexity, fa, hasConcave, m);
    } 

    // rotating the patch to have the first edge as the first element of the list
    patchRotationRightToEdge(patch, patchConvexity);

    // getting the number of edges in the patch
    int edge = 0;
    for (int convexity : patchConvexity){
        if (convexity >= 1)
            edge++;
    }
    return edge;
}
