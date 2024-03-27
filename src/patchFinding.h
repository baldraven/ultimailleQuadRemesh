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

inline int bfs(int startFacet, FacetAttribute<int>& facetAttributes, Quads& mesh){
    std::queue<int> facetQueue;
    int defectCount = 0;
    std::vector<int> defectVertices;

    facetQueue.push(startFacet);
    while (!facetQueue.empty() && defectCount < 3){
        Halfedge facetHalfedge = Facet(mesh, facetQueue.front()).halfedge().prev();
        facetQueue.pop();

        // explore each vertex of the facet        
        for (int i = 0; i < FACET_VERTICES; i++){
            facetHalfedge = facetHalfedge.next();

            if (facetHalfedge.from().on_boundary())
                return -1;

            Halfedge exploringHalfedge = facetHalfedge;

            if (isDefect(exploringHalfedge.from(), defectVertices)) {
                defectCount++;
            }

            // explore the facets surrounding a given vertex to mark them
            for (int i = 0; i < MAX_VALENCE; i++){
                facetAttributes[exploringHalfedge.facet()] = 1;
                
                Halfedge otherHalfedge = exploringHalfedge.opposite().next();
                if (otherHalfedge.active() && facetAttributes[otherHalfedge.facet()] == 0){
                    facetQueue.push(otherHalfedge.facet());
                    exploringHalfedge = otherHalfedge;
                }
            }
        }
    }
    return 1;
}

 inline int getPatch(Facet facet, FacetAttribute<int>& facetAttributes, std::list<int>& halfedgePatch, std::list<int>& patchConvexity){
    if (facetAttributes[facet] == 0){
        throw std::runtime_error("Error: facet not in patch!");
    }

    // find a halfedge in the border of the patch, inside
    Halfedge currentHalfedge = facet.halfedge();
    while(facetAttributes[currentHalfedge.facet()] >= 1){
        currentHalfedge = currentHalfedge.next().next().opposite();
    }
    currentHalfedge = currentHalfedge.opposite();

    halfedgePatch.clear();
    patchConvexity.clear();

    // cycle through the patch to get all the halfedges inside a double linked list
    halfedgePatch.push_back(currentHalfedge);
    Halfedge startHalfedge = currentHalfedge;
    bool firstIteration = true;

    while (currentHalfedge != startHalfedge || firstIteration){
        // boundary issue
        if (currentHalfedge.opposite().prev().opposite() == -1){
            return -1;
        }

        currentHalfedge = currentHalfedge.opposite();
        for (int i=0; i<MAX_VALENCE; i++){
            currentHalfedge = currentHalfedge.prev().opposite();
            if (facetAttributes[currentHalfedge.facet()] >= 1){
                halfedgePatch.push_back(currentHalfedge);
                patchConvexity.push_back(i-1);
                break;
            }
        }

        firstIteration = false;
    }  

    halfedgePatch.pop_front();
    return 1;
}


inline void addConcaveFaces(std::list<int>& patch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, bool& hasConcave, Quads& m){
    hasConcave = false;
    for (auto [a, b] : zip(patch, patchConvexity)) {
        if (b < 0){
            Facet f = Halfedge(m, a).opposite().facet();
            fa[f] = 2;
            hasConcave = true;
        }
    }
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

inline int completingPatch(FacetAttribute<int>& fa, Quads& m, Vertex v, std::list<int>& patch, std::list<int>& patchConvexity){
    if (v==15431)
        std::cout << "debug" << std::endl;
    bool hasConcave = true;
    while(hasConcave){
        int a = getPatch(v.halfedge().facet(), fa, patch, patchConvexity);
        if (a == -1)
            return -1;
        else if (a == 1000)
            return 1000;
        addConcaveFaces(patch, patchConvexity, fa, hasConcave, m);
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
