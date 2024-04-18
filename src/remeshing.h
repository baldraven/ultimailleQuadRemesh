#include <list>
#include <ultimaille/all.h>
#include <numeric>
#include <vector>
#include "ultimaille/algebra/vec.h"
#include "matrixEquations.h"
#include "ultimaille/surface.h"
#include "bvh.h"

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
    // TODO: put to true
    m.compact(false); 
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

inline int neoPuttingPointsInDefectEdges(int* partSegments, int iBarycentre, std::vector<int>& edgesAndDefectPointsOnBoundary, Quads& m, BVH bvh){
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

inline int getConstructedPointId(int i, int j, int* partSegments, Quads& m){
    int n = 0;
    for (int k=0; k<i-3; k=k+2){
        n += partSegments[k];
    }
    return m.nverts()-2-(n+j);
}


inline void meshingSubpatch(int* partSegments, int iBarycentre, std::list<int>& patch, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary, int& facet, BVH bvh){
    int n = edgesAndDefectPointsOnBoundary.size();

    // Construction the branches going from the barycentre to the points on the patch boundary
    int pointsAdded = neoPuttingPointsInDefectEdges(partSegments, iBarycentre, edgesAndDefectPointsOnBoundary, m, bvh);

    for (int i=1; i<n; i=i+2){
        
        int lines = partSegments[(i-2+n)%n];
        int columns = partSegments[i-1];
        std::cout << "i: " << i<<" columns:" << columns << " lines:" << lines << std::endl;
        
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
        //  - - - - -
        // | X X X X |
        // | X X X X |
        // | X X X X |
        // - - - - - -

            for (int j=1; j<columns; j++){
                // WIP: Intersection of 4 points instead of dividing in equal numbers
                
                // we create the new point and add it to the mesh
                // xi = x0 + i*(x1-x0)/n
                // x0 = point on the defect, x1 on the boundary

/*                 std::cout << "RUN: " << i << " " << k << " " << j << std::endl;
                std::cout << getConstructedPointId(i, k-1, partSegments, m) << std::endl; */

                
                vec3 x1 = Vertex(m, getConstructedPointId(i, k-1, partSegments, m)).pos(); // TO VALIDATE
                std::cout << "midpoint: i=" << i << " point=" << getConstructedPointId(i, k-1, partSegments, m) <<  std::endl;
                
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
        //  - - - - -
        // | + + + + X
        // | + + + + X
        // | + + + + X
        // - - - - -

            m.vert(facet, 1) = getConstructedPointId(i, k-1, partSegments, m);
            m.vert(facet-columns, 2) = getConstructedPointId(i, k-1, partSegments, m);
        }
    ///////////////////////////////////////////////////////

    //////////////////////////////////// Linking last line
        //  X X X X X
        // | + + + + |
        // | + + + + |
        // | + + + + |
        //  - - - - -

        it--;
        m.vert(facet-columns+1,3) = edgesAndDefectPointsOnBoundary[i-2];

        for (int j=1; j<columns; j++){
            int point = 0;
            if (i==1){
                point = getConstructedPointId(10, j-1, partSegments, m); // to validate
                std::cout << "I = 1: " << point<<std::endl;
            }
            else{
                point = getConstructedPointId((i-1)%10, j-1, partSegments, m);
                std::cout << "I > 1: " << point<<std::endl;
            }

            if (point==0){ // happens when only one facet line in the subpatch
                std::cout << "CASE POINT = 0" << std::endl;
                point = m.nverts()-pointsAdded-1;
            }

            m.vert(facet+j+1-columns, 3) = point; // TO VALIDATE
            m.vert(facet+j-columns, 2) = point;
        }
        m.vert(facet, 2) = iBarycentre;
        facet += columns;
    //////////////////////////////////////////
    }
}

inline void meshingSubpatchOld(int* partSegments, int iBarycentre, std::list<int>& patch, Quads& m, std::vector<int>& edgesAndDefectPointsOnBoundary, int& facet, BVH bvh){
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
        //  - - - - -
        // | X X X X |
        // | X X X X |
        // | X X X X |
        // - - - - - -

            for (int j=1; j<columns; j++){

                // WIP: Intersection of 4 points instead of dividing in equal numbers
    
                
                // we create the new point and add it to the mesh

                // xi = x0 + i*(x1-x0)/n
                // x0 = point on the defect, x1 on the boundary
                vec3 x1 = newPoints[k-1];
                vec3 x0 = boundaryPoint.pos();
                // We're gonna *cheat* and calculate the upper point of the defect edge, because it's unexistant for the first iteration.
                // TODO: redesign totally this function



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

inline bool remeshingPatch(std::list<int>& patch, std::list<int>& patchConvexity, int nEdge, Quads& m, FacetAttribute<int>& fa, int v, BVH bvh){
    assert(patchConvexity.front() >= 1);
    assert(nEdge == 3 || nEdge == 5);

    // DOES That actually works if she size is bigger than needed ?
    int segments[] = {0,0,0,0,0};
    int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
    segmentConstruction(patchConvexity, segments, nEdge);

    if (nEdge == 3){
        if (!solve3equations(segments, partSegments)){
            return false;
        }
    }
    else if (nEdge == 5){
        if (!solve5equations(segments, partSegments)){
            return false;
        }
    }

    std::cout << "solve" << nEdge << "equations success, root: " << v << std::endl;

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


    // TODO: change position or structure
    std::vector<int> edgesAndDefectPointsOnBoundary(2*nEdge,0);
    edgesAndDefectPointsOnBoundaryConstruction(patch, edgesAndDefectPointsOnBoundary, m, partSegments);

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

    // remove old facets and points (caution: it changes m topology)
    cleaningTopology(m, fa);

    return true;
}
