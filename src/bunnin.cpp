#include "helpers.h"
#include <cstdlib>
#include <iostream>
#include <ultimaille/all.h>
#include <list>

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

Facet getFacetById(int id, Quads& m){
    //efficiency ?
    for (auto f : m.iter_facets()){
        if (f == id)
            return f;
    }
    std::cout << "ERROR: facet not found !" << std::endl;
    exit(EXIT_FAILURE);
}

Halfedge getHalfedgeById(int id, Quads& m){
    //efficiency ?
    for (auto he : m.iter_halfedges()){
        if (he == id)
            return he;
    }
    std::cout << "ERROR: halfedge not found !" << std::endl;
    exit(EXIT_FAILURE);
}

int get_valence(Vertex v){
    int val = 0;
    Halfedge he = v.halfedge();
    Halfedge start = he;
    for (int i=0; i<8; i++){ //8 as maximum valence
        val++; 
        Halfedge otherHe = he.opposite().next();
        if (otherHe.active()){
            he = otherHe;
            if(he == start)
                return val;
        }
    }
    std::cout << "valence calculation error" << std::endl;
    exit(EXIT_FAILURE);
}

void bfs(int f, FacetAttribute<int>& fa, Quads& m){
    std::queue<int> q;
    int defectCount = 0;
    std::vector<int> defectsVerts;

    q.push(f);
    while (!q.empty()&&defectCount<3){
        Halfedge preloopHe = getFacetById(q.front(), m).halfedge().prev();
        q.pop();

        // explore each vertex of the facet        
        for (int i=0; i<4; i++){
            preloopHe = preloopHe.next();
            Halfedge he = preloopHe;

            // check if not in array of defects
            // get_valence might be redundant and thus inneficient
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
}

// might be trouble around boundaries
void getPatch(Facet f, FacetAttribute<int>& fa, std::list<int>& HePatch, std::list<int>& patchConvexity){
    if (fa[f] == 0){
        std::cout << "ERROR: facet not in patch !" << std::endl;
        return;
    }

    // find an halfedge in the border of the patch, inside
    Halfedge he = f.halfedge();
    while(fa[he.facet()] >= 1){
        he = he.next().next().opposite();
    }
    he = he.opposite();

    // cycle through the patch to get all the halfedges inside a double linked list
    HePatch.push_back(he);
    Halfedge heStart = he;
    bool firstIter = true;
    while (he != heStart || firstIter){
        //std::cout << "from: " << he.from() << " to: " << he.to() << " facet: " << he.facet() << std::endl;
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

bool addConcaveFaces(std::list<int>& HePatch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, bool& hasConcave, Quads& m){
    int patchSize = HePatch.size();
    hasConcave = false;
    for (int i=0; i<patchSize; i++){
        if (patchConvexity.front() < 0){
            Facet f = getHalfedgeById(HePatch.front(), m).opposite().facet();
            fa[f] = 2;
            hasConcave = true;

        }
        HePatch.pop_front();
        patchConvexity.pop_front();
    }
    return hasConcave;
}

// https://www.mcs.anl.gov/~fathom/meshkit-docs/html/Mesh_8cpp_source.html (Jaal)
int remesh5patch(const int *segments, int *partSegments){
    double M[10][10], rhs[10];

    //  Equations:
    //      b0 -a2   = 0
    //      b1 -a3   = 0
    //      b2 -a4   = 0
    //      b3 -a0   = 0
    //      b4 -a1   = 0
    //      a0 + b0  = s0
    //      a1 + b1  = s1
    //      a2 + b2  = s2
    //      a3 + b3  = s3
    //      a4 + b4  = s4

    // For more details; See Guy Bunin's paper.
    // M =
    //   a0  a1 a2  a3  a4  b0   b1 b2   b3  b4
    //   0   0  -1   0   0   1   0   0   0   0
    //   0   0   0  -1   0   0   1   0   0   0
    //   0   0   0   0  -1   0   0   1   0   0
    //  -1   0   0   0   0   0   0   0   1   0
    //   0  -1   0   0   0   0   0   0   0   1
    //   1   0   0   0   0   1   0   0   0   0
    //   0   1   0   0   0   0   1   0   0   0
    //   0   0   1   0   0   0   0   1   0   0
    //   0   0   0   1   0   0   0   0   1   0
    //   0   0   0   0   1   0   0   0   0   1

    // octave:38> inv(M)
    // ans =
    //  -0.5   0.5   0.5  -0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    //  -0.5  -0.5   0.5   0.5  -0.5   0.5   0.5  -0.5  -0.5   0.5
    //  -0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    //   0.5  -0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    //   0.5   0.5  -0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    //   0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    //   0.5   0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    //   0.5   0.5   0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    //  -0.5   0.5   0.5   0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    //  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5   0.5

    //   First Row
    //  -0.5   0.5   0.5  -0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    M[0][0] = -0.5;
    M[0][1] = 0.5;
    M[0][2] = 0.5;
    M[0][3] = -0.5;
    M[0][4] = -0.5;
    M[0][5] = 0.5;
    M[0][6] = -0.5;
    M[0][7] = -0.5;
    M[0][8] = 0.5;
    M[0][9] = 0.5;
    //  -0.5  -0.5   0.5   0.5  -0.5   0.5   0.5  -0.5  -0.5   0.5
    M[1][0] = -0.5;
    M[1][1] = -0.5;
    M[1][2] = 0.5;
    M[1][3] = 0.5;
    M[1][4] = -0.5;
    M[1][5] = 0.5;
    M[1][6] = 0.5;
    M[1][7] = -0.5;
    M[1][8] = -0.5;
    M[1][9] = 0.5;
    //  -0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    M[2][0] = -0.5;
    M[2][1] = -0.5;
    M[2][2] = -0.5;
    M[2][3] = 0.5;
    M[2][4] = 0.5;
    M[2][5] = 0.5;
    M[2][6] = 0.5;
    M[2][7] = 0.5;
    M[2][8] = -0.5;
    M[2][9] = -0.5;

    //   0.5  -0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    M[3][0] = 0.5;
    M[3][1] = -0.5;
    M[3][2] = -0.5;
    M[3][3] = -0.5;
    M[3][4] = 0.5;
    M[3][5] = -0.5;
    M[3][6] = 0.5;
    M[3][7] = 0.5;
    M[3][8] = 0.5;
    M[3][9] = -0.5;

    //   0.5   0.5  -0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    M[4][0] = 0.5;
    M[4][1] = 0.5;
    M[4][2] = -0.5;
    M[4][3] = -0.5;
    M[4][4] = -0.5;
    M[4][5] = -0.5;
    M[4][6] = -0.5;
    M[4][7] = 0.5;
    M[4][8] = 0.5;
    M[4][9] = 0.5;

    //   0.5  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5
    M[5][0] = 0.5;
    M[5][1] = -0.5;
    M[5][2] = -0.5;
    M[5][3] = 0.5;
    M[5][4] = 0.5;
    M[5][5] = 0.5;
    M[5][6] = 0.5;
    M[5][7] = 0.5;
    M[5][8] = -0.5;
    M[5][9] = -0.5;

    //   0.5   0.5  -0.5  -0.5   0.5  -0.5   0.5   0.5   0.5  -0.5
    M[6][0] = 0.5;
    M[6][1] = 0.5;
    M[6][2] = -0.5;
    M[6][3] = -0.5;
    M[6][4] = 0.5;
    M[6][5] = -0.5;
    M[6][6] = 0.5;
    M[6][7] = 0.5;
    M[6][8] = 0.5;
    M[6][9] = -0.5;

    //   0.5   0.5   0.5  -0.5  -0.5  -0.5  -0.5   0.5   0.5   0.5
    M[7][0] = 0.5;
    M[7][1] = 0.5;
    M[7][2] = 0.5;
    M[7][3] = -0.5;
    M[7][4] = -0.5;
    M[7][5] = -0.5;
    M[7][6] = -0.5;
    M[7][7] = 0.5;
    M[7][8] = 0.5;
    M[7][9] = 0.5;

    //  -0.5   0.5   0.5   0.5  -0.5   0.5  -0.5  -0.5   0.5   0.5
    M[8][0] = -0.5;
    M[8][1] = 0.5;
    M[8][2] = 0.5;
    M[8][3] = 0.5;
    M[8][4] = -0.5;
    M[8][5] = 0.5;
    M[8][6] = -0.5;
    M[8][7] = -0.5;
    M[8][8] = 0.5;
    M[8][9] = 0.5;

    //  -0.5  -0.5   0.5   0.5   0.5   0.5   0.5  -0.5  -0.5   0.5
    M[9][0] = -0.5;
    M[9][1] = -0.5;
    M[9][2] = 0.5;
    M[9][3] = 0.5;
    M[9][4] = 0.5;
    M[9][5] = 0.5;
    M[9][6] = 0.5;
    M[9][7] = -0.5;
    M[9][8] = -0.5;
    M[9][9] = 0.5;

    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = 0.0;
    rhs[3] = 0.0;
    rhs[4] = 0.0;

    rhs[5] = segments[0];
    rhs[6] = segments[1];
    rhs[7] = segments[2];
    rhs[8] = segments[3];
    rhs[9] = segments[4];

    // calcul de x la matrice résolue
    std::vector<int> x(10);
    for (int i = 0; i < 10; i++) {
        double sum = 0.0;
        for (int j = 0; j < 10; j++)
            sum += M[i][j] * rhs[j];
        x[i] = (int) sum;
    }

    std::cout << "1" << std::endl;
    for (int i = 0; i < 10; i++){
        std::cout << "solution: "<< x[i] << std::endl;
        if (x[i] < 1) return 0;
    }

    std::cout << "2" << std::endl;
    if (x[0] != x[8]) return 0;
    if (x[0] + x[5] != rhs[5]) return 0;
    partSegments[0] = x[0];
    partSegments[1] = x[5];

    std::cout << "3" << std::endl;
    if (x[1] != x[9]) return 0;
    if (x[1] + x[6] != rhs[6]) return 0;
    partSegments[2] = x[1];
    partSegments[3] = x[6];

    std::cout << "4" << std::endl;
    if (x[2] != x[5]) return 0;
    if (x[2] + x[7] != rhs[7]) return 0;
    partSegments[4] = x[2];
    partSegments[5] = x[7];

    std::cout << "5" << std::endl;
    if (x[3] != x[6]) return 0;
    if (x[3] + x[8] != rhs[8]) return 0;
    partSegments[6] = x[3];
    partSegments[7] = x[8];

    std::cout << "6" << std::endl;
    if (x[4] != x[7]) return 0;
    if (x[4] + x[9] != rhs[9]) return 0;
    partSegments[8] = x[4];
    partSegments[9] = x[9];

    return 1;
}

int main() {
    std::string path = getAssetPath();
    Quads m;
    //read_by_extension(path + "catorus_voxel_bound_smooth.geogram", m);
    read_by_extension(path + "cowhead.geogram", m);
    m.connect();

    DisjointSet ds(m.ncorners());

    for (auto f : m.iter_facets()){
        for (auto he : f.iter_halfedges()){
            if (!he.opposite().active())
                continue;
            ds.merge(he, he.opposite());
        }
    }

    // Get associate facet id to group id
    std::vector setIds = std::vector<int>(m.ncorners());
    ds.get_sets_id(setIds);

    // Get a facets
    int k = 0;
    FacetAttribute<int> fa(m, 0);
    for (auto v: m.iter_vertices()){
        FacetAttribute<int> fa(m, 0);
        if (get_valence(v) != 4){
            k++;
            //if k < 2, do nothing
            //k7 dosent work
      /*       if (k < 7){
                continue;
            } */
         
            bfs(v.halfedge().facet(), fa, m);
            std::list<int> patch;
            std::list<int> patchConvexity;
            getPatch(v.halfedge().facet(), fa, patch, patchConvexity);

            /* int patchSize = patch.size();
            std::cout << "patch size: " << patchSize << std::endl;
            //print patch
            for (int i=0; i<patchSize; i++){
                Halfedge he = getHalfedgeById(patch.front(), m);
                patch.pop_front();
                std::cout << "From | to | facet: " << he.from() << " | " << he.to() << " | " << he.facet() << std::endl;
            } */
            
            // iter through patch
            // if convexity <= 0, add face to fa
            // repeat until no more face to add
            // if convexity >= 1, mark as edge
            // return number of edges 
            bool hasConcave = false;
            addConcaveFaces(patch, patchConvexity, fa, hasConcave, m);
            while(hasConcave){
                getPatch(v.halfedge().facet(), fa, patch, patchConvexity);
                addConcaveFaces(patch, patchConvexity, fa, hasConcave, m);
            } 
            getPatch(v.halfedge().facet(), fa, patch, patchConvexity); // could do without that if addConcaveFaces wouldn't pop the list

            // edge is the number of halfedges in the patch that has convexity >= 1
            // we iterate the list patchConvexity to count the number of >=1 
            int edge = 0;
            int patchSize = patch.size();
            for (int i=0; i<patchSize; i++){
                if (patchConvexity.front() >= 1)
                    edge++;
                patchConvexity.pop_front();
            }
            // could be done cleaner 
            if (edge == 5){
                getPatch(v.halfedge().facet(), fa, patch, patchConvexity); // could do without that if addConcaveFaces wouldn't pop the list
          /*       for(int i=0; i<patchSize; i++){
                    std::cout << "patch convexity enum: " << patchConvexity.front() << std::endl;
                    patchConvexity.pop_front();
                } */
                int segments[] = {0,0,0,0,0};
                int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
                int edgeLength = 0;
                for (int i=0; i<patchSize; i++){
                    edgeLength++;
                    if (patchConvexity.front() >= 1){
                        segments[5-edge] = edgeLength;
                        edge--;
                        edgeLength = 0;
                    }
                    patchConvexity.pop_front();
                }
                std::cout << "segments: " << segments[0] << " " << segments[1] << " " << segments[2] << " " << segments[3] << " " << segments[4] << std::endl;

               // int segments[] = {4,2,3,3,6};

                if (remesh5patch(segments, partSegments)){
                    std::cout << "remesh5patch success" << std::endl;
                } else {
                    std::cout << "remesh5patch failed" << std::endl;
                }
            }
        
            if (k > 10 ){
                break;
            }
        }
    }

    write_by_extension("bunnin.geogram", m, {{}, {{"fa", fa.ptr}}, {}});
    return 0;
}