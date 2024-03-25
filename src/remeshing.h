#include <list>
#include <ultimaille/all.h>
#include <numeric>
#include "ultimaille/algebra/vec.h"
#include "matrixEquations.h"

using namespace UM;
using Halfedge = typename Surface::Halfedge;
using Facet = typename Surface::Facet;
using Vertex = typename Surface::Vertex;

////////////////////////////////////////////////////////////////////////////////////////////////////////
// Remeshing

inline vec3 getBarycentre(std::list<int>, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary){
    vec3 barycentre = vec3(0,0,0);
    for (int i=0; i<10; i++)
        barycentre += Vertex(m, edgesAndDefectPointsOnBoundary[i]).pos();
    return barycentre/10;
}

inline void edgesAndDefectPointsOnBoundaryConstruction(std::list<int>& patch, std::vector<int>& edgesAndDefectPointsOnBoundary, Quads& m, int* partSegments){
    edgesAndDefectPointsOnBoundary[0] = Halfedge(m, patch.front()).from();
    auto it = patch.begin();
    for (int i=1; i<10; i++){
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
    m.compact(true); 
}

inline void puttingPointsInDefectEdges(int nbPointsToDivide, int i, int iBarycentre, int& pointsAdded, std::vector<vec3>& newPoints, std::vector<int>& newPointsIndex, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary){
    for (int j=1; j < nbPointsToDivide; j++){
        vec3 x1 = m.points[iBarycentre];
        vec3 x0 = Vertex(m, edgesAndDefectPointsOnBoundary[i]).pos();
        vec3 newPoint = x0 + j*(x1-x0) / (nbPointsToDivide);
        newPoints[j-1] = newPoint;
        newPointsIndex[j-1] = m.nverts()-pointsAdded-1;
        m.points[m.nverts()-pointsAdded-1] = newPoint;
        pointsAdded++;
    }
}

inline void meshingSubpatch(int* partSegments, int iBarycentre, std::list<int>& patch, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary, int& facet){
    int pointsAdded = 1;
    int thatOneFacet = -1;

    for (int i=1; i<10; i=i+2){
        int lines = partSegments[(i-2+10)%10];
        int columns = partSegments[i-1];


        // putting points in the defect edges
        int nbPointsToDivide = partSegments[(i+1)%10];
        std::vector<vec3> newPoints(nbPointsToDivide, vec3{0,0,0});
        std::vector<int> newPointsIndex(nbPointsToDivide, 0);
        puttingPointsInDefectEdges(nbPointsToDivide, i, iBarycentre, pointsAdded, newPoints, newPointsIndex, m, edgesAndDefectPointsOnBoundary);

 
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
                m.points[newPointIndex] = newPoint;

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
        //  X X X X X
        // | + + + + |
        // | + + + + |
        // | + + + + |
        //  - - - - -

        if(i==9){
            m.vert(thatOneFacet-lines+1, 3) = edgesAndDefectPointsOnBoundary[9];
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

inline void remeshing5patch(std::list<int>& patch, std::list<int>& patchConvexity, Quads& m, FacetAttribute<int>& fa, int v){
    assert(patchConvexity.front() >= 1);

    int segments[] = {0,0,0,0,0};
    int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
    segmentConstruction(patchConvexity, segments, 5);

    if (!solve5equations(segments, partSegments)){
        std::cout << "solve5equations failed, root: " << v << std::endl; return;
        return;
    }

    std::cout << "solve5equations success, root: " << v << std::endl;
    
    // Ensure that segments are valid
    assert(segments[0] == partSegments[0] + partSegments[1]);
    assert(segments[1] == partSegments[2] + partSegments[3]);
    assert(segments[2] == partSegments[4] + partSegments[5]);
    assert(segments[3] == partSegments[6] + partSegments[7]);
    assert(segments[4] == partSegments[8] + partSegments[9]);

    std::vector<int> edgesAndDefectPointsOnBoundary(10,0);
    edgesAndDefectPointsOnBoundaryConstruction(patch, edgesAndDefectPointsOnBoundary, m, partSegments);

    // barycentre calculation for defect position
    
    // extra points and facets are getting cleaned up at the end
    m.disconnect();
    m.create_facets(500);
    m.connect();
    m.points.create_points(500);

    vec3 barycentre = getBarycentre(patch, m, edgesAndDefectPointsOnBoundary);
    int iBarycentre = m.nverts()-1;
    m.points[iBarycentre] = barycentre;

    int facet = m.nfacets()-499;

    // working each subset of the patch (the 5 quads delimited by the defect)
    meshingSubpatch(partSegments, iBarycentre, patch, m, edgesAndDefectPointsOnBoundary,facet);

    // remove old facets and points (caution: it changes m topology)
    cleaningTopology(m, fa);
}
