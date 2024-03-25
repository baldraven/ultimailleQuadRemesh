#include <list>
#include <ultimaille/all.h>

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// patch finding

inline int get_valence(Vertex v){
    int val = 0;
    Halfedge he = v.halfedge();
    Halfedge start = he;
    if (v==1804)
        std::cout << "start: " << he.from() << "|" << he.to() << std::endl;
    for (int i=0; i<8; i++){ //8 as maximum valence
        val++; 
        Halfedge otherHe = he.opposite().next(); // valgrind error here
      
        if (otherHe.active()){
            if (he.from() == he.to()){ // workaround for this specific valgrind error when iterating a second time over the maillage
                return 4;
            }

            he = otherHe;

            if(he == start){
                return val;
            }
        }
    }
    std::cout << "ERROR: valence calculation: "<< v << std::endl;
    exit(EXIT_FAILURE);
    return 0;
}

inline int bfs(int f, FacetAttribute<int>& fa, Quads& m){
    std::queue<int> q;
    int defectCount = 0;
    std::vector<int> defectsVerts;

    q.push(f);
    while (!q.empty()&&defectCount<3){
        Halfedge preloopHe = Facet(m, q.front()).halfedge().prev();
        q.pop();

        // explore each vertex of the facet        
        for (int i=0; i<4; i++){
            preloopHe = preloopHe.next();

            if (preloopHe.from().on_boundary())
                return -1;

            Halfedge he = preloopHe;

            if (get_valence(he.from()) != 4
            // seek into the array of defects if the vertex is already in
            && std::find(defectsVerts.begin(),defectsVerts.end(), he.from())==defectsVerts.end()){
                defectCount++;
                defectsVerts.push_back(he.from());
            }

            // explore the facets surrouding a given vertex to mark them
            for (int i=0; i<8; i++){ // 8 as maximum valence
                fa[he.facet()] = 1; // add facet in the attribute
                Halfedge otherHe = he.opposite().next();

                if (otherHe.active() && fa[otherHe.facet()] == 0){
                    q.push(otherHe.facet());
                    he = otherHe;
                }
            }
        }
    }
    return 1;
}

inline void getPatch(Facet f, FacetAttribute<int>& fa, std::list<int>& HePatch, std::list<int>& patchConvexity){
    if (fa[f] == 0){
        std::cout << "ERROR: facet not in patch !" << std::endl;
        exit(EXIT_FAILURE);
    }

    // find an halfedge in the border of the patch, inside
    Halfedge he = f.halfedge();
    while(fa[he.facet()] >= 1){
        he = he.next().next().opposite();
    }
    he = he.opposite();
    if (he.from().on_boundary()){
        he = he.next().opposite();
    }

    HePatch.clear();
    patchConvexity.clear();

    // cycle through the patch to get all the halfedges inside a double linked list
    HePatch.push_back(he);
    Halfedge heStart = he;
    bool firstIter = true;

    while (he != heStart || firstIter){
        he = he.opposite().prev().opposite();

        if (fa[he.facet()] >= 1){
            HePatch.push_back(he);
            patchConvexity.push_back(-1);
        } else if (fa[he.prev().opposite().facet()] >= 1){
            he = he.prev().opposite();
            HePatch.push_back(he);
            patchConvexity.push_back(0);
        } else if (fa[he.prev().opposite().prev().opposite().facet()] >= 1){
            he = he.prev().opposite().prev().opposite();
            HePatch.push_back(he);
            patchConvexity.push_back(1);
        } else if (fa[he.prev().opposite().prev().opposite().prev().opposite().facet()] >= 1){
            he = he.prev().opposite().prev().opposite().prev().opposite();
            HePatch.push_back(he);
            patchConvexity.push_back(2);
        } else if (fa[he.prev().opposite().prev().opposite().prev().opposite().prev().opposite().facet()] >= 1){
            he = he.prev().opposite().prev().opposite().prev().opposite().prev().opposite();
            HePatch.push_back(he);
            patchConvexity.push_back(3);
        } else if (fa[he.prev().opposite().prev().opposite().prev().opposite().prev().opposite().prev().opposite().facet()] >= 1){
            he = he.prev().opposite().prev().opposite().prev().opposite().prev().opposite().prev().opposite();
            HePatch.push_back(he);
            patchConvexity.push_back(4);
        } else if (fa[he.prev().opposite().prev().opposite().prev().opposite().prev().opposite().prev().opposite().prev().opposite().facet()] >= 1){
            he = he.prev().opposite().prev().opposite().prev().opposite().prev().opposite().prev().opposite().prev().opposite();
            HePatch.push_back(he);
            patchConvexity.push_back(5);
        } else {
            std::cout << "ERROR: unable to find patch !" << std::endl;
            exit(EXIT_FAILURE);
            return;
        }
        firstIter = false;
    }  
    HePatch.pop_front();
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
    bool done = false;
    while (!done){
        if (patchConvexity.front() < 1){
            patch.splice(patch.begin(), patch, std::prev(patch.end()));
            patchConvexity.splice(patchConvexity.begin(), patchConvexity, std::prev(patchConvexity.end()));
        } else {
            done = true;
        }
    }
}

inline int completingPatch(FacetAttribute<int>& fa, Quads& m, Vertex v, std::list<int>& patch, std::list<int>& patchConvexity){
    bool hasConcave = true;
    while(hasConcave){
        getPatch(v.halfedge().facet(), fa, patch, patchConvexity);
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
