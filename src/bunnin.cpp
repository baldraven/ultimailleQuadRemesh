#include "helpers.h"
#include "ultimaille/algebra/vec.h"
#include "ultimaille/attributes.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <ultimaille/all.h>
#include <list>

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

Facet getFacetById(int id, Quads& m){
    for (auto f : m.iter_facets()){
        if (f == id)
            return f;
    }
    std::cout << "ERROR: facet not found !" << std::endl;
    exit(EXIT_FAILURE);
}

Halfedge getHalfedgeById(int id, Quads& m){
    for (auto he : m.iter_halfedges()){
        if (he == id)
            return he;
    }
    std::cout << "ERROR: halfedge not found !" << std::endl;
    exit(EXIT_FAILURE);
}

Vertex getVertexById(int id, Quads& m){
    for (auto v : m.iter_vertices()){
        if (v == id)
            return v;
    }
    std::cout << "ERROR: vertex not found !" << std::endl;
    exit(EXIT_FAILURE);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

int get_valence(Vertex v){
    int val = 0;
    Halfedge he = v.halfedge();
    Halfedge start = he;
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

// TODO: check around the boundaries
void getPatch(Facet f, FacetAttribute<int>& fa, std::list<int>& HePatch, std::list<int>& patchConvexity){
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

void patchRotationToEdge(std::list<int>& patch, std::list<int>& patchConvexity){
    bool done = false;
    while (!done){
        if (patchConvexity.front() < 1){
            patch.splice( patch.end(), patch, patch.begin() );
            patchConvexity.splice( patchConvexity.end(), patchConvexity, patchConvexity.begin() );
        } else {
            done = true;
        }
    }
}

void patchRotationRightToEdge(std::list<int>& patch, std::list<int>& patchConvexity){
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


void addConcaveFaces(std::list<int>& patch, std::list<int>& patchConvexity, FacetAttribute<int>& fa, bool& hasConcave, Quads& m){
    hasConcave = false;
    for (auto [a, b] : zip(patch, patchConvexity)) {
        if (b < 0){
            Facet f = getHalfedgeById(a, m).opposite().facet();
            fa[f] = 2;
            hasConcave = true;
        }
    }
}

int getNbEdgesInPatch(std::list<int>& patchConvexity){
    int edge = 0;
    for (int convexity : patchConvexity){
        if (convexity >= 1)
            edge++;
    }
    return edge;
}

void segmentConstruction(std::list<int>& patchConvexity, int* segments, int edge){
    // construct array with the number of points between each edge
    int edgeLength = 0;
    int skip = 1;
    for (int convexity : patchConvexity){
        if (skip){
            skip = 0;
            continue;
        } 
        edgeLength++;
        if (convexity >= 1){
            segments[5-edge] = edgeLength;
            edge--;
            edgeLength = 0;
        }
    }
    segments[4] = edgeLength+1;
}

vec3 getBarycentre(std::list<int>, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary){
    vec3 barycentre = vec3(0,0,0);
    for (int i=0; i<10; i++)
        barycentre += getVertexById(edgesAndDefectPointsOnBoundary[i], m).pos();
    return barycentre/10;
}

void edgesAndDefectPointsOnBoundaryConstruction(std::list<int>& patch, std::vector<int>& edgesAndDefectPointsOnBoundary, Quads& m, int* partSegments){
    edgesAndDefectPointsOnBoundary[0] = getHalfedgeById(patch.front(), m).from();
    auto it = patch.begin();
    for (int i=1; i<10; i++){
        std::advance(it, partSegments[i-1]);
        edgesAndDefectPointsOnBoundary[i] = getHalfedgeById(*it, m).from();
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

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
    // inverse M
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

    M[0][0]=-0.5;M[0][1]=0.5;M[0][2]=0.5;M[0][3]=-0.5;M[0][4]=-0.5;M[0][5]=0.5;M[0][6]=-0.5;M[0][7]=-0.5;M[0][8]=0.5;M[0][9]=0.5;
    M[1][0]=-0.5;M[1][1]=-0.5;M[1][2]=0.5;M[1][3]=0.5;M[1][4]=-0.5;M[1][5]=0.5;M[1][6]=0.5;M[1][7]=-0.5;M[1][8]=-0.5;M[1][9]=0.5;
    M[2][0]=-0.5;M[2][1]=-0.5;M[2][2]=-0.5;M[2][3]=0.5;M[2][4]=0.5;M[2][5]=0.5;M[2][6]=0.5;M[2][7]=0.5;M[2][8]=-0.5;M[2][9]=-0.5;
    M[3][0]=0.5;M[3][1]=-0.5;M[3][2]=-0.5;M[3][3]=-0.5;M[3][4]=0.5;M[3][5]=-0.5;M[3][6]=0.5;M[3][7]=0.5;M[3][8]=0.5;M[3][9]=-0.5;
    M[4][0]=0.5;M[4][1]=0.5;M[4][2]=-0.5;M[4][3]=-0.5;M[4][4]=-0.5;M[4][5]=-0.5;M[4][6]=-0.5;M[4][7]=0.5;M[4][8]=0.5;M[4][9]=0.5;
    M[5][0]=0.5;M[5][1]=-0.5;M[5][2]=-0.5;M[5][3]=0.5;M[5][4]=0.5;M[5][5]=0.5;M[5][6]=0.5;M[5][7]=0.5;M[5][8]=-0.5;M[5][9]=-0.5;
    M[6][0]=0.5;M[6][1]=0.5;M[6][2]=-0.5;M[6][3]=-0.5;M[6][4]=0.5;M[6][5]=-0.5;M[6][6]=0.5;M[6][7]=0.5;M[6][8]=0.5;M[6][9]=-0.5;
    M[7][0]=0.5;M[7][1]=0.5;M[7][2]=0.5;M[7][3]=-0.5;M[7][4]=-0.5;M[7][5]=-0.5;M[7][6]=-0.5;M[7][7]=0.5;M[7][8]=0.5;M[7][9]=0.5;
    M[8][0]=-0.5;M[8][1]=0.5;M[8][2]=0.5;M[8][3]=0.5;M[8][4]=-0.5;M[8][5]=0.5;M[8][6]=-0.5;M[8][7]=-0.5;M[8][8]=0.5;M[8][9]=0.5;
    M[9][0]=-0.5;M[9][1]=-0.5;M[9][2]=0.5;M[9][3]=0.5;M[9][4]=0.5;M[9][5]=0.5;M[9][6]=0.5;M[9][7]=-0.5;M[9][8]=-0.5;M[9][9]=0.5;

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

    //std::cout << "Passing step 1" << std::endl;
    for (int i = 0; i < 10; i++){
        if (x[i] < 1) return 0;
    }

    //std::cout << "Passing step 2" << std::endl;
    if (x[0] != x[8]) return 0;
    if (x[0] + x[5] != rhs[5]) return 0;
    partSegments[0] = x[0];
    partSegments[1] = x[5];

    //std::cout << "Passing step 3" << std::endl;
    if (x[1] != x[9]) return 0;
    if (x[1] + x[6] != rhs[6]) return 0;
    partSegments[2] = x[1];
    partSegments[3] = x[6];

    //std::cout << "Passing step 4" << std::endl;
    if (x[2] != x[5]) return 0;
    if (x[2] + x[7] != rhs[7]) return 0;
    partSegments[4] = x[2];
    partSegments[5] = x[7];

    //std::cout << "Passing step 5" << std::endl;
    if (x[3] != x[6]) return 0;
    if (x[3] + x[8] != rhs[8]) return 0;
    partSegments[6] = x[3];
    partSegments[7] = x[8];

    //std::cout << "Passing step 6" << std::endl;
    if (x[4] != x[7]) return 0;
    if (x[4] + x[9] != rhs[9]) return 0;
    partSegments[8] = x[4];
    partSegments[9] = x[9];

    return 1;
}

void cleaningTopology(Quads& m, FacetAttribute<int>& fa){
    for (int i=0; i < m.nfacets(); i++){
        if (fa[i] > 0){
            m.conn.get()->active[i] = false;
        }
    }
    m.compact(true); 
}

void puttingPointsInDefectEdges(int nbPointsToDivide, int i, int iBarycentre, int& pointsAdded, std::vector<vec3>& newPoints, std::vector<int>& newPointsIndex, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary){
    for (int j=1; j < nbPointsToDivide; j++){
        vec3 x1 = m.points[iBarycentre];
        vec3 x0 = getVertexById(edgesAndDefectPointsOnBoundary[i], m).pos();
        vec3 newPoint = x0 + j*(x1-x0) / (nbPointsToDivide);
        newPoints[j-1] = newPoint;
        newPointsIndex[j-1] = m.nverts()-pointsAdded-1;
        m.points[m.nverts()-pointsAdded-1] = newPoint;
        pointsAdded++;
    }
}

void meshingSubpatch(int* partSegments, int iBarycentre, std::list<int>& patch, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary, int& pointsAdded, int& facet, int& thatOneFacet){

    for (int i=1; i<10; i=i+2){

        // putting points in the defect edges
        int nbPointsToDivide = partSegments[(i+1)%10];
        std::vector<vec3> newPoints(nbPointsToDivide, vec3{0,0,0});
        std::vector<int> newPointsIndex(nbPointsToDivide, 0);
        puttingPointsInDefectEdges(  // see how to reduce number of arguments
            nbPointsToDivide, i, iBarycentre, pointsAdded, newPoints,
            newPointsIndex, m, edgesAndDefectPointsOnBoundary
        );

        int lines = partSegments[(i-2+10)%10];
        int columns = partSegments[i-1];

        // we need to go to the first halfedge of the subpatch
        auto it = patch.begin();
        std::advance(it, std::accumulate(partSegments, partSegments + i - 1, 0)); 


        // line 1 is different because boundary so we do it first
        m.vert(facet+1, 0) = edgesAndDefectPointsOnBoundary[i-1];
        for (int j=1; j<columns; j++){
            it++;
            int botBoundaryPointIndex = getHalfedgeById(*it, m).from();
            m.vert(facet+j+1, 0) = botBoundaryPointIndex;
            m.vert(facet+j, 1) = botBoundaryPointIndex;
        }
        it++;
        facet += columns;
        m.vert(facet, 1) = getHalfedgeById(*it, m).from();

        //////////////////////////////////////////////////////////

        // remaining lines

        // we go back to the first halfedge of the subpatch
        std::advance(it, -columns -1 + (i!=1?1:0));

        for (int k=1; k<lines; k++){ 
            // first point is in the boundary so we do it apart
            it--;

            Vertex boundaryPoint = getHalfedgeById(*it, m).from();
            m.vert(facet+1, 0) = boundaryPoint;
            m.vert(facet+1-columns, 3) = boundaryPoint;
            /////////////////////////////

            for (int j=1; j<columns; j++){
                
                // we create the new point and add it to the mesh

                // xi = x0 + i*(x1-x0)/n
                // x0 = point on the defect, x1 on the boundary
                vec3 x1 = newPoints[k-1];
                vec3 x0 = boundaryPoint.pos();

                vec3 newPoint = x0 + j*(x1-x0)/columns;

                int newPointIndex = m.nverts()-pointsAdded-1;
                m.points[newPointIndex] = newPoint;

                // we link it to facets
                m.vert(facet+j+1, 0) = newPointIndex;
                m.vert(facet+j, 1) = newPointIndex;
                m.vert(facet-columns+j, 2) = newPointIndex;
                m.vert(facet-columns+j+1, 3) = newPointIndex;

                pointsAdded++;
            }
            facet += columns;
                


            // Linking last column to the boundary
            m.vert(facet, 1) = newPointsIndex[k-1];
            m.vert(facet-columns, 2) = newPointsIndex[k-1]; 
            ////////////////////////////

            // we need memory of the this facet because we want to link those facets in the last subpatch
            // because the points on the defect edge are created in the last subpatch, and it was supposed to be linked in the first one
            if (k == lines-1 && i == 1){
                thatOneFacet = facet;
            }

        }
        // we were doing it above but the loop might not get entered at all
        if (lines == 1 && i == 1){
                thatOneFacet = facet;
         }
        ////////////////// Linking last line
        if (i == 1) // the one of the first subpatch will be linked in the last subpatch
            continue; 

        it--;
        m.vert(facet-columns+1,3) = edgesAndDefectPointsOnBoundary[i-2];

        int prevLines = partSegments[(i-4+10)%10];
        int prevColumns = partSegments[i-3];
        int numberOfPointsInPrevSubpatch = (prevLines-1)*(prevColumns-1);

        for (int j=1; j<columns; j++){
            if (newPointsIndex[0] == 0) // happens when only one facet line in the subpatch
                newPointsIndex[0] = m.nverts()-pointsAdded-1;
            int boundaryPoint = newPointsIndex[0]+columns-j + numberOfPointsInPrevSubpatch;
            m.vert(facet+j+1-columns, 3) = boundaryPoint; 
            m.vert(facet+j-columns, 2) = boundaryPoint;
        }
        m.vert(facet, 2) = iBarycentre;
        facet += columns;
        //////////////////////////////////////////


        // Linking the problematic line, the one that should have been linked in the first subpatch
        if(i==9){
            m.vert(thatOneFacet-lines+1, 3) = edgesAndDefectPointsOnBoundary[9];
            for (int j=1; j<lines; j++){
                int boundaryPoint = newPointsIndex[lines-j-1];
                m.vert(thatOneFacet-j, 2) = boundaryPoint;
                m.vert(thatOneFacet-j+1, 3) = boundaryPoint;
            }
            m.vert(thatOneFacet, 2) = iBarycentre;
        } 
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////

int main() {
    std::string path = getAssetPath();
    Quads m;
    //read_by_extension(path + "catorus_voxel_bound_smooth.geogram", m);
    read_by_extension(path + "cowhead.geogram", m);
    m.connect();


    FacetAttribute<int> fa(m, 0);
    // iterating through the vertices until finding a defect 
    for (int o=0; o < 4; o++){
        std::cout << "-----------------o: " << o << std::endl;
    for (auto v: m.iter_vertices()){

        if (v>248)
            break;  
        
        if (get_valence(v) == 4)
            continue;
        
        // reset the attribute
        fa.fill(0);
        
        // constructing the patch and coloring the facets inside
        bfs(v.halfedge().facet(), fa, m);

        // expanding the patch to include concave facets. patch is a list of halfedges in the boundary of the patch
        std::list<int> patch;
        std::list<int> patchConvexity;
        bool hasConcave = true;
        while(hasConcave){
            getPatch(v.halfedge().facet(), fa, patch, patchConvexity);
            addConcaveFaces(patch, patchConvexity, fa, hasConcave, m);
        } 
        
        // computing the number of edges in the patch
        int edge = getNbEdgesInPatch(patchConvexity);

        // Currently we only work with pentagons
        if (edge == 5){

            // rotating the patch to have the first edge as the first element of the list
            patchRotationRightToEdge(patch, patchConvexity);
            assert(patchConvexity.front() >= 1);

            int segments[] = {0,0,0,0,0};
            int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
            segmentConstruction(patchConvexity, segments, edge);

            if (!remesh5patch(segments, partSegments)){
                std::cout << "remesh5patch failed, root: " << v << std::endl;
                continue;
            }
            std::cout << "remesh5patch success, root: " << v << std::endl;
            
            //std::cout << "partSegments: " << partSegments[0] << " " << partSegments[1] << " " << partSegments[2] << " " << partSegments[3] << " " << partSegments[4] << " " << partSegments[5] << " " << partSegments[6] << " " << partSegments[7] << " " << partSegments[8] << " " << partSegments[9] << std::endl;
            if (v==22)
                continue;
           
            // Ensure that segments are valid
            assert(segments[0] == partSegments[0] + partSegments[1]);
            assert(segments[1] == partSegments[2] + partSegments[3]);
            assert(segments[2] == partSegments[4] + partSegments[5]);
            assert(segments[3] == partSegments[6] + partSegments[7]);
            assert(segments[4] == partSegments[8] + partSegments[9]);

            std::vector<int> edgesAndDefectPointsOnBoundary(10,0);
            edgesAndDefectPointsOnBoundaryConstruction(patch, edgesAndDefectPointsOnBoundary, m, partSegments);

            // barycentre calculation for defect position
            vec3 barycentre = getBarycentre(patch, m, edgesAndDefectPointsOnBoundary);
            
            // extra points and facets are getting cleaned up at the end
            m.disconnect();
            m.create_facets(500);
            m.connect();

            m.points.create_points(400);
            int iBarycentre = m.nverts()-1;
            m.points[iBarycentre] = barycentre;

            int pointsAdded = 1;
            int facet = m.nfacets()-499;
            int thatOneFacet = -1;
          

            // working each subset of the patch (the 5 quads delimited by the defect)
            meshingSubpatch(
                partSegments, iBarycentre, patch, m, edgesAndDefectPointsOnBoundary,
                pointsAdded, facet, thatOneFacet
            );

            // remove old facets and points (caution: it changes m topology)
            cleaningTopology(m, fa);
        }
        }
    }

    write_by_extension("bunnin.geogram", m, {{}, {{"patch", fa.ptr}}, {}});
    return 0;
}
