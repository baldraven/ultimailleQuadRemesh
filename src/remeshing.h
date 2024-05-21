#include <cmath>
#include <list>
#include <ostream>
#include <ultimaille/all.h>
#include <numeric>
#include <unistd.h>
#include <vector>
#include "ultimaille/algebra/vec.h"
#include "matrixEquations.h"
#include "ultimaille/attributes.h"
#include "ultimaille/surface.h"
#include "bvh.h"
#include <assert.h>


using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remeshing

inline vec3 getBarycentre(std::list<int>, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary){
    int n = edgesAndDefectPointsOnBoundary.size();
    vec3 barycentre = vec3(0,0,0);
    for (int i=0; i<n; i++)
        barycentre += Vertex(m, edgesAndDefectPointsOnBoundary[i]).pos();
    return barycentre/n;
}

inline void edgesAndDefectPointsOnBoundaryConstruction(std::list<int>& patch, std::vector<int>& edgesAndDefectPointsOnBoundary, Quads& m, int* partSegments){
    int n = edgesAndDefectPointsOnBoundary.size();
    edgesAndDefectPointsOnBoundary[0] = Halfedge(m, patch.front()).from();
    auto it = patch.begin();
    for (int i=1; i<n; i++){
        std::advance(it, partSegments[i-1]);
        edgesAndDefectPointsOnBoundary[i] = Halfedge(m, *it).from();
    }
}

inline void cleaningTopology(Quads& m, FacetAttribute<int>& fa){
    for (int i=0; i < m.nfacets(); i++){
        if (fa[i] > 0){
            m.conn.get()->active[i] = false;
        }
    }
    //std::cout << "NOCOMPACT" << std::endl;
    
    m.compact(true); 
}

inline void puttingPointsInDefectEdges(int nbPointsToDivide, int i, int iBarycentre, int& pointsAdded, std::vector<vec3>& newPoints, std::vector<int>& newPointsIndex, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary, BVH bvh){
    for (int j=1; j < nbPointsToDivide; j++){
        vec3 x1 = m.points[iBarycentre];
        vec3 x0 = Vertex(m, edgesAndDefectPointsOnBoundary[i]).pos();
        vec3 newPoint = x0 + j*(x1-x0) / (nbPointsToDivide);
        newPoints[j-1] = newPoint;
        newPointsIndex[j-1] = m.nverts()-pointsAdded-1;
        m.points[m.nverts()-pointsAdded-1] = bvh.project(newPoint);
        pointsAdded++;
    }
}

inline int neoPuttingPointsInDefectEdges(int* partSegments, int iBarycentre, std::vector<int>& edgesAndDefectPointsOnBoundary, Quads& m){
    int n = edgesAndDefectPointsOnBoundary.size();
    int pointsAdded = 2; //only the barycentre for now
    int arraySize = 0;
    for (int i=0; i<10; i=i+2){
        arraySize += partSegments[i];
    }
    arraySize -= n/2;
    std::cout << "ArraySize: "<<arraySize << std::endl;
    

    for (int i=1; i<10; i+=2){
        std::cout << "Initial I: " << i << std::endl;
        
        int nbPointsToDivide = partSegments[(i+1)%n];

        for (int j=1; j < nbPointsToDivide; j++){
            vec3 x1 = m.points[iBarycentre];
            vec3 x0 = Vertex(m, edgesAndDefectPointsOnBoundary[i]).pos();
            vec3 newPoint = x0 + j*(x1-x0) / (nbPointsToDivide);
            std::cout << "THE POINT: " << m.nverts()-pointsAdded<<std::endl;
            
            m.points[m.nverts()-pointsAdded] = newPoint;
            pointsAdded++;
        }
    }

    assert(pointsAdded == arraySize+2);
    return pointsAdded;
}


inline void divideInSubPatches(int* partSegments, int iBarycentre, std::vector<int>& edgesAndDefectPointsOnBoundary, Quads& m, std::list<int>& patch){
    int n = edgesAndDefectPointsOnBoundary.size();

    int arraySize = 0;
    for (int i=0; i<10; i=i+2){
        arraySize += partSegments[i];
    }
    arraySize -= n/2;
    m.points.create_points(arraySize);

    int pointsAdded = 0;

    for (int i=1; i<10; i+=2){
        
        int nbPointsToDivide = partSegments[(i+1)%n];

        for (int j=1; j < nbPointsToDivide; j++){
            vec3 x1 = m.points[iBarycentre];
            vec3 x0 = Vertex(m, edgesAndDefectPointsOnBoundary[i]).pos();
            vec3 newPoint = x0 + j*(x1-x0) / (nbPointsToDivide);
            
            m.points[m.nverts()-pointsAdded] = newPoint;
            pointsAdded++;
        }
    }


    std::vector<std::list<int>> subpatches(n/2, std::list<int>());
    for (int i=0; i<n/2; i++){
        // edgesAndDefectPointsOnBoundary : e
        // ei -> ei+1
        // ei+1 -> bi
        // b(i-1mod10) -> e(i-1mod10)
        // e(i-1mod10) -> ei
        for (int j=0; j<partSegments[i*2]; j++){
            //subpatches[i].push_back();
        }

    }


}


inline Triangles quand2tri(Quads& m){
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

inline void meshingQuad(std::vector<int>& anodes, std::vector<int>& bnodes, std::vector<int>& cnodes, std::vector<int>& dnodes, Quads& m, BVH bvh){
    // assert that we have a topological rectangle
    assert(anodes.size() == cnodes.size());
    assert(bnodes.size() == dnodes.size());

    // assert that edges of the rectangle are correctly connected and oriented
    //        b ->
    //      - - - - -
    //   ^ | + + + + | ^
    // a | | + + + + | | c
    //     | + + + + |  
    //     - - - - - -
    //       d ->
    assert(anodes[anodes.size()-1] == bnodes[0]);
    assert(bnodes[bnodes.size()-1] == cnodes[cnodes.size()-1]);
    assert(dnodes[dnodes.size()-1] == cnodes[0]);
    assert(anodes[0] == dnodes[0]);

    int a = anodes.size();
    int b = bnodes.size();

    m.points.create_points((a-2)*(b-2));

    for (int i=1; i<a; i++){
        for (int j=1; j<b; j++){
            //  - - - - -
            // | X X X X |
            // | X X X X |
            // | X X X X |
            // - - - - - -

            // we'll construct the points starting from the bottom left and decrementing starting from m.nverts()
            // TODO: consider the 4 points on the grid instead of just 2
            // TODO: projecting

            vec3 x0 = Vertex(m, anodes[i]).pos();
            vec3 x1 = Vertex(m, cnodes[i]).pos();

            vec3 newPoint = x0 + j*(x1-x0)/ (b-1); 

            int newPointIndex = m.nverts()-  ((i-1)*(b-2) + j);
            int btmNewPointIndex = newPointIndex + (b-2); 
            if (i<a-1 && j<b-1)
                m.points[newPointIndex] = bvh.project(newPoint);


            if (i==1 && j==1)
                m.conn->create_facet({anodes[1], newPointIndex, dnodes[1], dnodes[0]});
            else if (i==1 && j<b-1)
                m.conn->create_facet({newPointIndex+1, newPointIndex, dnodes[j], dnodes[j-1]});
            else if (i==1 && j==b-1)
                m.conn->create_facet({newPointIndex+1, cnodes[i], dnodes[j], dnodes[j-1]});
            else if (i<a-1 && j==1)
                m.conn->create_facet({anodes[i], newPointIndex, btmNewPointIndex, anodes[i-1]});
            else if (i<a-1 && j<b-1)
                m.conn->create_facet({newPointIndex+1, newPointIndex, btmNewPointIndex, btmNewPointIndex+1});
            else if (i<a-1 && j==b-1)
                m.conn->create_facet({newPointIndex+1, cnodes[i], cnodes[i-1], btmNewPointIndex+1});
            else if (i==a-1 && j==1)
                m.conn->create_facet({bnodes[0], bnodes[1], btmNewPointIndex, anodes[a-2]});
            else if (i==a-1 && j<b-1)
                m.conn->create_facet({bnodes[j-1], bnodes[j], btmNewPointIndex, btmNewPointIndex+1});
            else if (i==a-1 && j==b-1)
                m.conn->create_facet({bnodes[b-2], cnodes[a-1], cnodes[a-2], btmNewPointIndex+1});
        }
    }
}

inline void meshingSubpatch(int* partSegments, int iBarycentre, std::list<int>& patch, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary, int& facet, BVH bvh){
    int pointsAdded = 1;
    int thatOneFacet = -1;

    int n = edgesAndDefectPointsOnBoundary.size();

    // Ensure that segments are valid

    for (int i=1; i<n; i=i+2){
        int lines = partSegments[(i-2+n)%n];
        int columns = partSegments[i-1];


        // putting points in the defect edges
        int nbPointsToDivide = partSegments[(i+1)%n];
        std::vector<vec3> newPoints(nbPointsToDivide, vec3{0,0,0});
        std::vector<int> newPointsIndex(nbPointsToDivide, 0);
        puttingPointsInDefectEdges(nbPointsToDivide, i, iBarycentre, pointsAdded, newPoints, newPointsIndex, m, edgesAndDefectPointsOnBoundary, bvh);

 
        // we need to go to the first halfedge of the subpatch
        auto it = patch.begin();
        std::advance(it, std::accumulate(partSegments, partSegments + i - 1, 0)); 


        // Connecting and creating all points in the subpatch

/////////// line 1 is different because boundary so we do it first
        // X are the points being connected
        //
        //  - - - - -
        // | + + + + |
        // | + + + + |
        // | + + + + |
        // X X X X X X

        m.vert(facet+1, 0) = edgesAndDefectPointsOnBoundary[i-1];
        for (int j=1; j<columns; j++){
            it++;
            int botBoundaryPointIndex = Halfedge(m, *it).from();
            m.vert(facet+j+1, 0) = botBoundaryPointIndex;
            m.vert(facet+j, 1) = botBoundaryPointIndex;
        }
        it++;
        facet += columns;
        m.vert(facet, 1) = Halfedge(m, *it).from();

////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////// remaining lines

        // we go back to the first halfedge of the subpatch
        std::advance(it, -columns -1 + (i!=1?1:0));

        for (int k=1; k<lines; k++){ 

    //// first points are on the boundary so we do it apart
        // X are the points being connected
        //
        //  - - - - -
        // X + + + + |
        // X + + + + |
        // X + + + + |
        // - - - - - -

            it--;
            Vertex boundaryPoint = Halfedge(m, *it).from();
            m.vert(facet+1, 0) = boundaryPoint;
            m.vert(facet+1-columns, 3) = boundaryPoint;


    ///////////////////////////////////////////////////////

    /////////////////////////////////// points in the middle
        // X are the points being connected
        //
        //  - - - - -
        // | X X X X |
        // | X X X X |
        // | X X X X |
        // - - - - - -

            for (int j=1; j<columns; j++){

    
                
                // we create the new point and add it to the mesh

                // xi = x0 + i*(x1-x0)/n
                // x0 = point on the defect, x1 on the boundary
                vec3 x1 = newPoints[k-1];
                vec3 x0 = boundaryPoint.pos();

                vec3 newPoint = x0 + j*(x1-x0)/columns;

                int newPointIndex = m.nverts()-pointsAdded-1;
                m.points[newPointIndex] = bvh.project(newPoint);

                // we link it to facets
                m.vert(facet+j+1, 0) = newPointIndex;
                m.vert(facet+j, 1) = newPointIndex;
                m.vert(facet-columns+j, 2) = newPointIndex;
                m.vert(facet-columns+j+1, 3) = newPointIndex;

                pointsAdded++;
            }
            facet += columns;
    ///////////////////////////////////////////////////////         

    /////////////////////////////////////////// Last column
        // X are the points being connected
        //
        //  - - - - -
        // | + + + + X
        // | + + + + X
        // | + + + + X
        // - - - - -

            m.vert(facet, 1) = newPointsIndex[k-1];
            m.vert(facet-columns, 2) = newPointsIndex[k-1]; 

    ///////////////////////////////////////////////////////

            // we need memory of the this facet because we want to link those facets in the last subpatch
            // because the points on the defect edge are created in the last subpatch, and it was supposed to be linked in the first one
            if (k == lines-1 && i == 1){
                thatOneFacet = facet;
            }

        }
        // we need memory of the this facet because we want to link those facets in the last subpatch
        // we were doing it above but the loop might not get entered at all
        if (lines == 1 && i == 1){
                thatOneFacet = facet;
        }

    //////////////////////////////////// Linking last line
        // X are the points being connected
        //
        //  X X X X X
        // | + + + + |
        // | + + + + |
        // | + + + + |
        //  - - - - -


        // the one of the first subpatch will be linked in the last subpatch
        if (i == 1) 
            continue; 

        it--;
        m.vert(facet-columns+1,3) = edgesAndDefectPointsOnBoundary[i-2];

        int prevLines = partSegments[(i-4+n)%n];
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
        //  X X X X X
        // | + + + + |
        // | + + + + |
        // | + + + + |
        //  - - - - -

        if(i==n-1){
            m.vert(thatOneFacet-lines+1, 3) = edgesAndDefectPointsOnBoundary[n-1];
            for (int j=1; j<lines; j++){
                int boundaryPoint = newPointsIndex[lines-j-1];
                m.vert(thatOneFacet-j, 2) = boundaryPoint;
                m.vert(thatOneFacet-j+1, 3) = boundaryPoint;
            }
            m.vert(thatOneFacet, 2) = iBarycentre;
        } 
    //////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
    }
}

inline void segmentConstruction(std::list<int>& patchConvexity, int* segments, int edge){
    // construct array with the number of points between each edge
    int n = edge;
    int edgeLength = 0;
    int skip = 1;
    for (int convexity : patchConvexity){
        if (skip){
            skip = 0;
            continue;
        } 
        edgeLength++;
        if (convexity >= 1){
            segments[n-edge] = edgeLength;
            edge--;
            edgeLength = 0;
        }
    }
    segments[n-1] = edgeLength+1;
}

inline void fillingConvexPos(std::list<int>& patchConvexity, std::vector<int>& convexPos){
    int count = 0;
    int letterToFill = 0;
    for (int v : patchConvexity){
        if (v >= 1){
            convexPos[letterToFill]=count;
            letterToFill++;
        }
        count++;
    }
}

inline int pyMod(int a, int b){
    return (a%b+b)%b;
}

inline bool testRotations(std::vector <int>& convexPos, std::vector<int>& cumulConvexity, int& rotation, int size){
    rotation = 0;
    for (int i=0; i<(int)convexPos.size();i++){
        rotation = convexPos[i];
        if (pyMod(convexPos[i]-rotation,size)==cumulConvexity[0]
            && pyMod(convexPos[(i+1)%4]-rotation,size)==cumulConvexity[1]
            && pyMod(convexPos[(i+2)%4]-rotation,size)==cumulConvexity[2]
            && pyMod(convexPos[(i+3)%4]-rotation,size)==cumulConvexity[3])
        {
            return true;
        }
    }
    return false;
}

inline void rotateToA(std::list<int>& patch, std::list<int>& patchConvexity, int a, int b, int c, int d){
    std::vector<int> convexPos = {0,0,0,0};
    std::vector<int> cumulConvexity = {0,0,0,0};
    cumulConvexity[1] = a;
    cumulConvexity[2] = a+b;
    cumulConvexity[3] = a+b+c; 

    // Filling convex pos
    fillingConvexPos(patchConvexity, convexPos);
    
    int rotation = 0;
    int size = patch.size();
    bool found = testRotations(convexPos, cumulConvexity, rotation, size);

    if (!found){
        patch.reverse();
        patchConvexity.reverse();
        fillingConvexPos(patchConvexity, convexPos);
        found = testRotations(convexPos, cumulConvexity, rotation, size);
    }
    assert(found);

    // rotating the patch
    for (int i=0; i<rotation; i++){
        patch.push_back(patch.front());
        patch.pop_front();
        patchConvexity.push_back(patchConvexity.front());
        patchConvexity.pop_front();
    }
}


inline int roundUpDivide(int a, int b){
    return (a+b-1)/b;
}

inline bool find(std::list<int>& v, int x){
    return std::find(v.begin(), v.end(), x) != v.end();
}

inline bool solve4equations(int* segments, std::list<int>& patch, std::list<int>& patchConvexity, Quads& m, FacetAttribute<int>& fa){
    if (segments[0] == segments[2] && segments[1] == segments[3]){
        std::cout << "PERFECT QUAD REMESH POSSIBLE" << std::endl;
        return true;
    }
    int a = fmax(segments[0], segments[2]);
    int c = fmin(segments[0], segments[2]);
    int b = fmin(segments[1], segments[3]);
    int d = fmax(segments[1], segments[3]);

    int segmentsTri[] = {d-b,  c, a};
    int partSegments[] = {0,0,0,0,0,0};
    if (solve3equations(segmentsTri, partSegments)){
        //std::cout << "SOLUTION 1 FOUND !" << std::endl;

        // First we'll rotate the patch so the patch starts at beginning of a followed by b
        rotateToA(patch, patchConvexity, a, b, c, d);
        write_by_extension("outputA.geogram", m, {{}, {{"patch", fa.ptr}, }, {}});
        
        // TODO: add .from()'s 

        a++;b++;c++;d++; // because we want to include the last points
        // Then we'll construct the left regular quad
        std::vector<int> anodes(a, 0);
        auto it = patch.begin();
        for (int i=0; i<a; i++){
            anodes[i] = Halfedge(m, *it);
            it++;
        }

        // b nodes goes from the last node of a (included) to half of b
        std::vector<int> bnodes(roundUpDivide(b, 2), 0);
        it--;
        for (int i=0; i<roundUpDivide(b, 2); i++){
            bnodes[i] = Halfedge(m, *it);
            it++;
        }

        // c nodes must be constructed
        std::vector<int> cnodes(a, 0);
        it --;
        Halfedge h = Halfedge(m, *it);
        cnodes[0] = h.from();

        it = patch.end();
        std::advance(it, -b/2);
        cnodes[a-1] = Halfedge(m, *it).from();

        // we construct a-2 points distributed between cnodes[0] and cnodes[a-1] and put them in cnodes
        m.points.create_points(a-2);
        for (int i=1; i<a-1; i++){
            // TODO: verify point positions
            vec3 x0 = Vertex(m, cnodes[0]).pos();
            vec3 x1 = Vertex(m, cnodes[a-1]).pos();
            vec3 newPoint = x0 + i*(x1-x0)/(a-1);
            cnodes[i] = m.nverts()-i;
            m.points[m.nverts()-i] = newPoint;
        }
        std::reverse(cnodes.begin(), cnodes.end());


        // Let's do d nodes now
        std::vector<int> dnodes(roundUpDivide(b, 2), 0);
        it = patch.begin();
        it++;
        for (int i=0; i<roundUpDivide(b, 2); i++){
            dnodes[i] = Halfedge(m, *it).from();
            it--;
        }

        exit(0);

        return false;
    }




    // TODO: Sideway triangle insertion
    //int segmentsTri2[] = {d, b, a-c};
    //if (solve3equations(segmentsTri2, partSegments)){
    //    std::cout << "SOLUTION 2 FOUND !" << std::endl;
    //    return false;
    //}


    // Remesh left part 
    // Remesh right part
    // Remesh triangle :)
    // il va nous falloir un *patch* pour chacune des parties
    // algorithme: on parcourt le patch jusqu'à atteindre le point (dans le cas 1) en d-b -c ou qq chose comme ça
    // on fait un sharp turn pour avoir un angle de 1 vers l'intérieur
    // on compte une nouvelle fois le nombre de points (c) et on fait un sharp turn pour avoir un angle de 1 vers l'intérieur
    // on referme

    return false;
}

inline void patchToNodes(std::list<int>& patch, int* segments, std::vector<int>& anodes, std::vector<int>& bnodes, std::vector<int>& cnodes, std::vector<int>& dnodes, Quads& m){
    int aSize = segments[1]+1;
    int bSize = segments[0]+1;

    auto it = patch.begin();   
    for (int i = 0; i < (int) patch.size(); i++){

        if (i==0){
            anodes[aSize-1] = Halfedge(m, *it).from();
            bnodes[0] = Halfedge(m, *it).from();
        } else if (i < segments[0]){
            bnodes[i] = Halfedge(m, *it).from();
        } else if (i == segments[0]){
            assert(i==bSize-1);
            bnodes[bSize-1] = Halfedge(m, *it).from();
            cnodes[aSize-1] = Halfedge(m, *it).from();
        } else if (i < segments[0]+segments[1]){
            cnodes[aSize-1-i+segments[0]] = Halfedge(m, *it).from();
        } else if (i == segments[0]+segments[1]){
            cnodes[0] = Halfedge(m, *it).from();
            dnodes[bSize-1] = Halfedge(m, *it).from();
        } else if (i < segments[0]+segments[1]+segments[2]){
            dnodes[bSize-1-i+segments[0]+segments[1]] = Halfedge(m, *it).from();
        } else if (i == segments[0]+segments[1]+segments[2]){
            dnodes[0] = Halfedge(m, *it).from();
            anodes[0] = Halfedge(m, *it).from();
        } else {
            anodes[i-segments[0]-segments[1]-segments[2]] = Halfedge(m, *it).from();
        }
        it++;
    }
    std::cout << "anodes: ";
    for (int i : anodes){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << "bnodes: ";
    for (int i : bnodes){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << "cnodes: ";
    for (int i : cnodes){
        std::cout << i << " ";
    }
    std::cout << std::endl;
    std::cout << "dnodes: ";
    for (int i : dnodes){
        std::cout << i << " ";
    }
    std::cout << std::endl;
}

inline bool remeshingPatch(std::list<int>& patch, std::list<int>& patchConvexity, int nEdge, Quads& m, FacetAttribute<int>& fa, int v, BVH bvh){
    assert(patchConvexity.front() >= 1);
    assert(nEdge == 3 || nEdge == 5 || nEdge == 4);

    int segments[] = {0,0,0,0,0};
    int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
    segmentConstruction(patchConvexity, segments, nEdge);


    if (nEdge == 4){
        //print segments
        //std::cout << "Segments: " << segments[0] << " " << segments[1] << " " << segments[2] << " " << segments[3] << std::endl;
        if (!solve4equations(segments, patch, patchConvexity, m, fa)){
            return false;
        }
    }
    else if (nEdge == 3){
        if (!solve3equations(segments, partSegments)){
            return false;
        }
    }
    else if (nEdge == 5){
        if (!solve5equations(segments, partSegments)){
            return false;
        }
    }

    // patch into triangles for projection
    // first we're gonna have to create a new mesh with just the patch
    // TODO: Put that in a function
    Quads mPatch;
    std::vector<int> facetsInPatch;
    for (int i = 0; i < m.nfacets(); i++){
        if (fa[i] > 0){
            facetsInPatch.push_back(i);
        }
    }
    // Deep copying mesh 
    mPatch.points.create_points(m.nverts());
    for (Vertex v : m.iter_vertices()){
        mPatch.points[v] = v.pos();
    }
    mPatch.create_facets(facetsInPatch.size());
    for (int i = 0; i < (int) facetsInPatch.size(); i++){
        Facet f = Facet(m, facetsInPatch[i]);
        mPatch.vert(i, 0) = m.vert(f, 0);
        mPatch.vert(i, 1) = m.vert(f, 1);
        mPatch.vert(i, 2) = m.vert(f, 2);
        mPatch.vert(i, 3) = m.vert(f, 3);
    }
    m.compact(true);
    Triangles mTri = quand2tri(mPatch);
    BVH bvhPatch(mTri);


    std::cout << "solve" << nEdge << "equations success, root: " << v << std::endl;

    std::vector<int> edgesAndDefectPointsOnBoundary(2*nEdge,0);
    edgesAndDefectPointsOnBoundaryConstruction(patch, edgesAndDefectPointsOnBoundary, m, partSegments);


    


    if (nEdge == 3 || nEdge == 5){
        // extra points and facets are getting cleaned up at the end
        // TODO: Optimize facet number creation
        m.disconnect();
        m.create_facets(2000);
        int facet = m.nfacets()-1999;
        m.connect();
        m.points.create_points(2000);


        // barycentre calculation for defect position
        vec3 barycentre = getBarycentre(patch, m, edgesAndDefectPointsOnBoundary);
        int iBarycentre = m.nverts()-1;
        m.points[iBarycentre] = bvhPatch.project(barycentre);

        // working each subset of the patch (the 5 quads delimited by the defect)
        meshingSubpatch(partSegments, iBarycentre, patch, m, edgesAndDefectPointsOnBoundary,facet, bvhPatch);




    } else if (nEdge == 4){
        
        for (auto i : patch){
            std::cout << Halfedge(m, i).from() << " ";
        }
        std::cout << std::endl;


        // Constructiong the sides
        int aSize = segments[1]+1;
        int bSize = segments[0]+1;
        std::vector<int> anodes(aSize);
        std::vector<int> bnodes(bSize);
        std::vector<int> cnodes(aSize);
        std::vector<int> dnodes(bSize);
        
        patchToNodes(patch, segments, anodes, bnodes, cnodes, dnodes, m);

        meshingQuad(anodes, bnodes, cnodes, dnodes, m, bvhPatch);
    }
    


    // remove old facets and points (caution: it changes m topology)
    cleaningTopology(m, fa);
  /*   if (nEdge==4){
        std::cout << "Vs: ";
        write_by_extension("outputQuad.geogram", m);
        exit(0);
    } */


    return true;
}
