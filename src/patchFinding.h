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

    for (int i = 0; i < 6; i++){
        valence++; 
        Halfedge nextHalfedge = currentHalfedge.opposite().next();
      
        if (nextHalfedge.active()){
            currentHalfedge = nextHalfedge;

            if (currentHalfedge == startHalfedge)
                return valence;
        }
    }

    return 4;
/*  throw std::runtime_error("Error in valence calculation for vertex: " + std::to_string(vertex));
    return 0; */
}

inline bool isDefect(Vertex v, std::vector<int>& defects) {
    if (getValence(v) != 4 && 
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
    Halfedge facetHalfedge = Facet(mesh, startFacet).halfedge();
    Halfedge exploringHalfedge = facetHalfedge;
    Halfedge otherHalfedge = exploringHalfedge;

    while (!facetQueue.empty() && defectCount < 3){
        facetHalfedge = Facet(mesh, facetQueue.front()).halfedge().prev();
        facetQueue.pop();

        // explore each vertex of the facet        
        for (int i = 0; i < 4; i++){
            facetHalfedge = facetHalfedge.next();

            if (facetHalfedge.from().on_boundary())
                return -1;

            exploringHalfedge = facetHalfedge;

            if (isDefect(exploringHalfedge.from(), defectVertices)){
                defectCount++;

                // marking the facets surrounding the last defect vertex
                if (defectCount == 3){
                    for (int i = 0; i < MAX_VALENCE; i++){
                        facetAttributes[exploringHalfedge.facet()] = 1;
                        exploringHalfedge = exploringHalfedge.opposite().next();
                    }
                    break;
                }
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

    // find a halfedge in the border of the patch, inside the patch
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
    int nbIteration = 0;

    int MAX_ITERATION = 500;

    while ((boundaryHe != startHalfedge || firstIteration) && nbIteration < MAX_ITERATION){
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

        nbIteration++;
        firstIteration = false;
    }  

    halfedgePatch.pop_front();
    return 1;
}

inline int updateBoundaryHe(int& boundaryHe, Halfedge& he, Quads& m, FacetAttribute<int>& fa){
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

// TODO: Put patch updating in a function ?
inline int expandPatch(int& boundaryHe, std::list<int>& patch, FacetAttribute<int>& fa, Quads& m){
    Halfedge he = Halfedge(m, boundaryHe);
    for (int i : patch) {
        he = Halfedge(m, i).opposite();

        if (fa[he.next().next().opposite().facet()] >= 1){
            return -1;
        }

        fa[he.facet()] = 2;
    }

    return updateBoundaryHe(boundaryHe, he, m, fa);
}


inline int addConcaveFaces(int& boundaryHe, std::list<int>& patch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, bool& hasConcave, Quads& m){
    hasConcave = false;
    Halfedge he = Halfedge(m, boundaryHe);
    for (auto [a, b] : zip(patch, patchConvexity)) {
        if (b < 0){           

            he = Halfedge(m, a).opposite();

            if (fa[he.next().next().opposite().facet()] >= 1){
            //std::cout << "This is reached" << std::endl;
                return -1;
            }

            if (fa[he.facet()] >= 1){
               //std::cout << "This is reached" << std::endl;
                return -1;
            }


            fa[he.facet()] = 2;

            hasConcave = true;
        }
    }

    return updateBoundaryHe(boundaryHe, he, m, fa);
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

inline void animateDebug2(Quads& m, int i, FacetAttribute<int>& fa, int iter){
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
    std::string s = "../animation/outputDEBUG" + std::to_string(iter) + number + ".geogram";
    write_by_extension(s, m, {{}, {{"patch", fa.ptr}, }, {}});
}

inline int completingPatch(int boundaryHe, FacetAttribute<int>& fa, Quads& m, std::list<int>& patch, std::list<int>& patchConvexity, int t, int v, int iter){
    // Unmarking the facets of the patch
    for (int i : patch){
        if (fa[Halfedge(m, i).facet()] == 3)
            fa[Halfedge(m, i).facet()] = 2;
    }

    bool hasConcave = true;
    for (int i : patchConvexity){
        if (i >= 1){
            hasConcave = false;
            break;
        }
    }

    // TODO: What if it was the first attempt and the patch was not concave ? in this case we don't want to expand
    if (!hasConcave){
        if(expandPatch(boundaryHe, patch, fa, m) == -1)
            return -1;
        hasConcave = true;
    }

    while(hasConcave ){
        if (getPatch(Halfedge(m, boundaryHe), fa, patch, patchConvexity) == -1)
            return -1;


        if (addConcaveFaces(boundaryHe, patch, patchConvexity, fa, hasConcave, m) == -1)
            return -1;

    } 


    // Mark the outline of the patch
    for (int i : patch){
        fa[Halfedge(m, i).facet()] = 3;
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

    // TODO: we still have a topological disk problem sometimes (v=191)
}
