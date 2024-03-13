#include "helpers.h"
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <ultimaille/all.h>
#include <list>

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;



// template usage ?
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

Vertex getVertexById(int id, Quads& m){
    //efficiency ?
    for (auto v : m.iter_vertices()){
        if (v == id)
            return v;
    }
    std::cout << "ERROR: vertex not found !" << std::endl;
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
            if(he == start){
                return val;
            }
        }
    }
    std::cout << "ERROR: valence calculation" << std::endl;
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

    std::cout << "Passing step 1" << std::endl;
    for (int i = 0; i < 10; i++){
        if (x[i] < 1) return 0;
    }

    std::cout << "Passing step 2" << std::endl;
    if (x[0] != x[8]) return 0;
    if (x[0] + x[5] != rhs[5]) return 0;
    partSegments[0] = x[0];
    partSegments[1] = x[5];

    std::cout << "Passing step 3" << std::endl;
    if (x[1] != x[9]) return 0;
    if (x[1] + x[6] != rhs[6]) return 0;
    partSegments[2] = x[1];
    partSegments[3] = x[6];

    std::cout << "Passing step 4" << std::endl;
    if (x[2] != x[5]) return 0;
    if (x[2] + x[7] != rhs[7]) return 0;
    partSegments[4] = x[2];
    partSegments[5] = x[7];

    std::cout << "Passing step 5" << std::endl;
    if (x[3] != x[6]) return 0;
    if (x[3] + x[8] != rhs[8]) return 0;
    partSegments[6] = x[3];
    partSegments[7] = x[8];

    std::cout << "Passing step 6" << std::endl;
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


    Quads m2;





    FacetAttribute<int> fa(m, 0);

    // iterating through the vertices until finding a defect 
    for (auto v: m.iter_vertices()){

        if (get_valence(v) == 4)
            continue;

        fa.fill(0);
        // valence calculation error near k=150 to debug
        
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
        int edge = 0;
        for (int convexity : patchConvexity){
            if (convexity >= 1)
                edge++;
        }

        // Currently we only work with pentagons
        if (edge == 5){

            // constructing sides length data for remeshing
            int segments[] = {0,0,0,0,0};
            int edgeLength = 0;
            for (int convexity : patchConvexity){
                edgeLength++;
                if (convexity >= 1){
                    segments[5-edge] = edgeLength;
                    edge--;
                    edgeLength = 0;
                }
            }

            int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
            if (remesh5patch(segments, partSegments)){
                std::cout << "remesh5patch success" << std::endl;
        
                // print partsegment
                std::cout << "partSegments: " << partSegments[0] << " " << partSegments[1] << " " << partSegments[2] << " " << partSegments[3] << " " << partSegments[4] << " " << partSegments[5] << " " << partSegments[6] << " " << partSegments[7] << " " << partSegments[8] << " " << partSegments[9] << std::endl;


                // getting the points that will be attached to the defect
                std::vector<int> simpleRemeshPoints(10,0);
                simpleRemeshPoints[0] = getHalfedgeById(patch.back(), m).from();
                auto it = patch.begin();
                it--;
                for (int i=1; i<10; i++){
                    for (int j=0; j<partSegments[i-1]; j++){
                        it++;
                    }
                    simpleRemeshPoints[i] = getHalfedgeById(*it, m).from();
                }


                // barycentre calculation for defect position
                vec3 barycentre = vec3(0,0,0);
                for (int i=0; i<10; i++){
                    barycentre += getVertexById(simpleRemeshPoints[i], m).pos();
                }
                barycentre /= 10;
                std::cout << "barycentre: " << barycentre << std::endl;


                // constructing the new mesh with single defect
                m.points.create_points(100);
                int iBarycentre = m.nverts()-1;
                m.points[iBarycentre] = barycentre;
                m.create_facets(5);

                // linking the unique defect to the mesh -- we don't have to do that if we do the submeshing
                int firstRemeshFacet = m.nfacets()-5;
                m.vert(firstRemeshFacet, 0) = simpleRemeshPoints[0];
                m.vert(firstRemeshFacet, 1) = simpleRemeshPoints[1];
                m.vert(firstRemeshFacet, 2) = m.nverts()-1;
                m.vert(firstRemeshFacet, 3) = simpleRemeshPoints[9];
                int j = 1;
                for (int i=1; i<5; i++){
                    m.vert(firstRemeshFacet+i, 0) = simpleRemeshPoints[j];
                    m.vert(firstRemeshFacet+i, 1) = simpleRemeshPoints[j+1];
                    m.vert(firstRemeshFacet+i, 2) = simpleRemeshPoints[j+2];
                    m.vert(firstRemeshFacet+i, 3) = m.nverts()-1;
                    j=j+2;
                }

                // linking the other points
                int nbPoints = 0;
                // TODO: ajuster le nombre de points 
                nbPoints = std::accumulate(partSegments, partSegments+10, nbPoints);
                nbPoints -= 6;
                m.points.create_points(nbPoints);


                std::vector<int> nbFacetsPerSubsection = {0,0,0,0,0}; // testedOK
                for (int i=0; i<10; i=i+2){
                    nbFacetsPerSubsection[i/2] = partSegments[(i-1+10)%10]*partSegments[i];
                }
                int sumFacets = std::accumulate(nbFacetsPerSubsection.begin(), nbFacetsPerSubsection.end(), 0);
                m.create_facets(500);

                int pointsAdded = 0;
                int facet = m.nfacets()+100-sumFacets; // should we add +1
                for (int i=1; i<10; i=i+2){
                    // TODO: impove readability
                    // putting points in the defect edges
                    int nbPointsToDivide = partSegments[(i+1)%10];
                    std::vector<vec3> newPoints(nbPointsToDivide, vec3{0,0,0});
                    for (int j=1; j < nbPointsToDivide; j++){
                        //xi = x0 + i*(x1-x0)/n
                        vec3 newPoint = m.points[iBarycentre] + j*(getVertexById(simpleRemeshPoints[i], m).pos() - m.points[iBarycentre]) / (nbPointsToDivide+1);
                        newPoints[j-1] = newPoint;
                        m.points[m.nverts()-pointsAdded-1] = newPoint;
                        pointsAdded++;
                    }
//############################################################################################################
                    // working in this patch subsection
                    // let's work line by line
                    int lines = partSegments[(i-1+10)%10];
                    int columns = partSegments[i];
                    auto it = patch.begin();

                    // line 1 is different becaues boundary so we do it first
                    //std::cout << "pointIt: " << simpleRemeshPoints[0] << std::endl;
                    m.vert(facet, 0) = simpleRemeshPoints[0];
                    for (int j=1; j<columns; j++){ //end bound to revisit
                        //std::cout << "pointIt: " << getHalfedgeById(*it, m).from() << std::endl;
                        m.vert(facet+j, 1) = getHalfedgeById(*it, m).from();
                        m.vert(facet+j+1, 0) = getHalfedgeById(*it, m).from();
                        it++;
                    }
                    facet += columns;

                     std::cout << "are newPoints actually legit: " << newPoints[0] << "   " << newPoints[1] << std::endl;
                    // remaining lines
                    it = patch.begin();
                    for (int k=1; k<lines; k++){
                        // first point is in the boundary so we do it apart
                        it--;
                        Vertex boundaryPoint = getHalfedgeById(*it, m).from();
                        m.vert(facet-columns, 3) = boundaryPoint;
                        m.vert(facet, 0) = boundaryPoint;
                        for (int j=1; j<columns; j++){
                            // create the new point (zzzzz)
                            //xi = x0 + i*(x1-x0)/n
                            //vec3 newPoint = m.points[iBarycentre] + j*(getVertexById(simpleRemeshPoints[i], m).pos() - m.points[iBarycentre]) / (nbPointsToDivide+1);
                            //m.points[m.nverts()-pointsAdded-1] = newPoint;
                            // x0 = point on the defect, x1 on the boundary
                            
                            // the points are not correct !!!!!!!!!

                            vec3 x0 = newPoints[j-1];
                            vec3 x1 = boundaryPoint.pos();
                            vec3 newPoint = x0 + k*(x1-x0)/lines;
                            std::cout << "x1 and x2: " << x0 << " " << x1 << " " << newPoint << std::endl;
                            m.points[m.nverts()-pointsAdded-1] = newPoint;
                            std::cout << "newpoint: " << m.nverts()-pointsAdded-1 << " " << newPoint << std::endl;
                            pointsAdded++;
                            // m.vert(facet-columns+j, 2) = 
                        }
                    }
                    break;
                }



                m.connect();

                // remove old facets (caution: it changes m topology)
                for (int i=0; i < m.nfacets(); i++){
                    if (fa[i] > 0){
                        m.conn.get()->active[i] = false;
                    }
                }
                m.compact(false); 

                break;
            } else {
                std::cout << "remesh5patch failed" << std::endl;
            }
        }
    }

    write_by_extension("bunnin.geogram", m, {{}, {{"patch", fa.ptr}}, {}});
    return 0;
}
