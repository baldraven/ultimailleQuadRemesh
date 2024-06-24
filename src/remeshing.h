#include <algorithm>
#include <iostream>
#include <iterator>
#include <list>
#include <ostream>
#include <ultimaille/all.h>
#include <utility>
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

inline int pyMod(int a, int b){
    return (a%b+b)%b;
}

inline void cleaningTopology(Quads& m, FacetAttribute<int>& fa){
    for (int i=0; i < m.nfacets(); i++){
        if (fa[i] > 0){
            m.conn.get()->active[i] = false;
        }
    }
    m.compact(true); 
}

inline void meshingRectangle(std::vector<int>& anodes, std::vector<int>& bnodes, std::vector<int>& cnodes, std::vector<int>& dnodes, Quads& m, BVH bvh){
    // TODO: add reverse as a parameter

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

    // In the case there's no new point to create, just have to connect the boundary of the rectangle
    if (anodes.size() < 3 && bnodes.size() < 3){
        m.conn->create_facet({anodes[0], bnodes[0], bnodes[1], cnodes[0]});
        return;
    }

    // If the rectangle is too small, we make the small part be on the columns rather than lines so we just have to deal with 1 case
    bool reversed = false;
    if (anodes.size()<3){
        std::swap(anodes, dnodes);
        std::swap(cnodes, bnodes);
        reversed = true;
    }
    int a = anodes.size();
    int b = bnodes.size();


    // Creating the new points inside the patch, starting from the bottob left corner
    m.points.create_points((a-2)*(b-2));
    for (int i=1; i<a; i++){
        for (int j=1; j<b; j++){

            int newPointIndex = 0;
            int btmNewPointIndex = 0;
            if (b<3){
                newPointIndex = cnodes[i];
                btmNewPointIndex = cnodes[i-1];
            }
            else{
                vec3 x0 = Vertex(m, anodes[i]).pos();
                vec3 x1 = Vertex(m, cnodes[i]).pos();

                vec3 newPoint = x0 + j*(x1-x0)/ (b-1); 

                newPointIndex = m.nverts()-  ((i-1)*(b-2) + j);
                btmNewPointIndex = newPointIndex + (b-2); 
                if (i<a-1 && j<b-1)
                    m.points[newPointIndex] = bvh.project(newPoint);
            }

            // Creating the facets with the new points. facets have an orientation so have a case where we reverse the order of the nodes to adjust
            if (reversed){
                if (i==1 && j==1)
                    m.conn->create_facet({dnodes[0], dnodes[1], newPointIndex, anodes[1]});
                else if (i==1 && j<b-1)
                    m.conn->create_facet({dnodes[j-1], dnodes[j], newPointIndex, newPointIndex+1});
                else if (i==1 && j==b-1)
                    m.conn->create_facet({dnodes[j-1], dnodes[j], cnodes[i], newPointIndex+1});
                else if (i<a-1 && j==1)
                    m.conn->create_facet({anodes[i-1], btmNewPointIndex, newPointIndex, anodes[i]});
                else if (i<a-1 && j<b-1)
                    m.conn->create_facet({btmNewPointIndex+1, btmNewPointIndex, newPointIndex, newPointIndex+1});
                else if (i<a-1 && j==b-1)
                    m.conn->create_facet({btmNewPointIndex+1, cnodes[i-1], cnodes[i], newPointIndex+1});
                else if (i==a-1 && j==1)
                    m.conn->create_facet({anodes[a-2], btmNewPointIndex, bnodes[1], bnodes[0]});
                else if (i==a-1 && j<b-1)
                    m.conn->create_facet({btmNewPointIndex+1, btmNewPointIndex, bnodes[j], bnodes[j-1]});
                else if (i==a-1 && j==b-1)
                    m.conn->create_facet({btmNewPointIndex+1, cnodes[a-2], cnodes[a-1], bnodes[b-2]});
            } else {
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
    if (reversed == true){
        std::swap(anodes, dnodes);
        std::swap(cnodes, bnodes);
    }
}

inline void constructBarycentre(int size, std::vector<std::vector<int>>& anodesList, Quads& m, BVH bvh, vec3& barycentrePos, int& barycentreIndex){
   std::vector<int> baryNodes(size, 0);
    for (int i=0;i<size;i++){
        baryNodes[i]=anodesList[i][anodesList[i].size()-1];
        barycentrePos += Vertex(m, baryNodes[i]).pos();
    }

    barycentrePos /= size;
    barycentrePos = bvh.project(barycentrePos);

    m.points.create_points(1);
    barycentreIndex = m.nverts()-1;
    m.points[barycentreIndex]=barycentrePos;
}

inline void nPatchRemesh(int* partSegments, std::list<int>& patch, Quads& m, int size, BVH bvh){
    // We have a 3 or a 5 patch, that we'll divide it in 3 or 5 rectangles to remesh them individually
    // each ones will have anodes, bnodes, cnodes, dnodes (see the meshingRectangle function)


    // anodes and dnodes, they are on the boundary on the initial patch so we iterate over it
    std::vector<std::vector <int>> anodesList(size);
    std::vector<std::vector <int>> dnodesList(size);
    auto it = patch.begin();
    it++;
    for (int i=0;i<size;i++){
        it--;
        for (int j=0;j<partSegments[2*i]+1;j++){
            anodesList[i].push_back(Halfedge(m, *it).from());
            it++;
        }
        it--;
        for (int j=0;j<partSegments[2*i+1]+1;j++){
            dnodesList[i].push_back(Halfedge(m, *it).from());
            it++;
        }
    }
    for (auto& nodeList : dnodesList){
        std::reverse(nodeList.begin(), nodeList.end());
    }
    std::rotate(dnodesList.rbegin(), dnodesList.rbegin() + 1, dnodesList.rend());
    dnodesList[0][0] = anodesList[0][0];


    int barycentreIndex = 0;
    vec3 barycentrePos = {0,0,0};
    constructBarycentre(size, anodesList, m, bvh, barycentrePos, barycentreIndex);

        
    // b nodes, we're going to create new points between the barycentre and the anodes
    std::vector<std::vector <int>> bnodesList(size);

    for (int i=0;i<size;i++){

        int x0index = anodesList[i][anodesList[i].size()-1];
        bnodesList[i].push_back(x0index);
        vec3 x0 = Vertex(m, x0index).pos();
        vec3 x1 = barycentrePos;
        int bnodesLen = (int)dnodesList[i].size();

        m.points.create_points(bnodesLen-1);
        for (int j=1;j<bnodesLen-1;j++){
            // Make the new point
            vec3 newPoint = x0 +j*(x1-x0) / (bnodesLen-1);
            int newPointIndex = m.nverts()-j;
            m.points[newPointIndex] = bvh.project(newPoint);
            bnodesList[i].push_back(newPointIndex);
        }
        bnodesList[i].push_back(barycentreIndex);
    } 

    // c nodes are the same as the bnodes of previous patch so we just rotate the list
    std::vector<std::vector <int>> cnodesList = bnodesList;
    std::rotate(cnodesList.rbegin(), cnodesList.rbegin() + 1, cnodesList.rend());


    for (int i=0; i<size; i++){
        meshingRectangle(anodesList[i], bnodesList[i], cnodesList[i], dnodesList[i], m, bvh);
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

inline bool rotateToA(std::list<int>& patch, std::list<int>& patchConvexity, int a, int b, int c){
    std::vector<int> convexPos = {0,0,0,0};
    std::vector<int> cumulConvexity = {0,0,0,0};
    cumulConvexity[1] = a;
    cumulConvexity[2] = a+b;
    cumulConvexity[3] = a+b+c; 

    fillingConvexPos(patchConvexity, convexPos);
    
    int rotation = 0;
    int size = patch.size();
    bool found = testRotations(convexPos, cumulConvexity, rotation, size);
    bool wasReversed = false;

    if (!found){
        patch.reverse();
        patchConvexity.reverse();
        fillingConvexPos(patchConvexity, convexPos);
        found = testRotations(convexPos, cumulConvexity, rotation, size);
        wasReversed = true;
    }
    assert(found);

    // rotating the patch
    for (int i=0; i<rotation; i++){
        patch.push_back(patch.front());
        patch.pop_front();
        patchConvexity.push_back(patchConvexity.front());
        patchConvexity.pop_front();
    }
    
    return wasReversed;
}

inline int roundUpDivide(int a, int b){
    return (a+b-1)/b;
}

inline bool find(std::list<int>& v, int x){
    return std::find(v.begin(), v.end(), x) != v.end();
}

inline void rectanglePatchRemesh(std::list<int>& patch, int* segments, Quads& m, BVH bvh){
    int aSize = segments[1]+1;
    int bSize = segments[0]+1;

    std::vector<int> anodes(aSize);
    std::vector<int> bnodes(bSize);
    std::vector<int> cnodes(aSize);
    std::vector<int> dnodes(bSize);

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
    meshingRectangle(anodes, bnodes, cnodes, dnodes, m, bvh);

}

inline void createPointsBetween2Vx(std::vector<int>& nodes, int n, Quads& m, BVH bvh){
    m.points.create_points(n-1);
    for (int i=1; i<n; i++){
        vec3 x0 = Vertex(m, nodes[0]).pos();
        vec3 x1 = Vertex(m, nodes[n]).pos();
        vec3 newPoint = x0 + i*(x1-x0)/n;
        nodes[i] = m.nverts()-i;
        m.points[m.nverts()-i] = bvh.project(newPoint);
    }
}

inline void ajustPartSegments(int* partSegments, int c, int btm, int a){
    std::vector<int> ints = {c, btm, a};
    for (int j=0; j<3 ; j++){
        bool changed = false;
        for (int i=0; i<3 ; i++){
            int side = partSegments[2*i]+partSegments[2*i+1];
            if (side != ints[i]){
                changed = true;
                std::rotate(partSegments, partSegments+2, partSegments+6);
                break;
            }
        }
        if (!changed){
            return;
        }
    }
    ints = {c, a, btm};
    for (int j=0; j<3 ; j++){
        bool changed = false;
        for (int i=0; i<3 ; i++){
            int side = partSegments[2*i]+partSegments[2*i+1];
            if (side != ints[i]){
                changed = true;
                std::rotate(partSegments, partSegments+2, partSegments+6);
                break;
            }
        }
        if (!changed){
            return;
        }
    }
    assert(false);
}

inline void quadrilateralPatchRemesh(int* partSegments, std::list<int>& patch, std::list<int>& patchConvexity, Quads& m, BVH bvh, int a, int b, int c, int d){
    //                 bnodes       bnodes2
    //               -------->    -------->
    //             ------------------------
    //           ^ |xxxxxxxxxxx|xxxxxxxxx|   ^
    //           | |xxxxxxxxxx|x|xxxxxxxx|   |
    // anodes    | |xxxxxxxx|x x x|xxxxxx|   | cnodes2
    //           | |xxxxxxx|x x x x|xxxxx|   |
    //           | |xxxxxx|x x x x x|xxxx|   |
    //             ------------------------
    //              -------->   ---->  ------->
    //                 dnodes   btmPart  dnodes2
    //
    //                   ||
    //        ^       | x x |         ^
    //        |      | x x x x |      |   anodes2
    // cnodes |    | x x x x x |      |
    //        |   | x x x x x x |     |
    //                ------------>
    //                   btmPart

    // First we'll rotate the patch so it starts at the bottom left corner
    bool wasReversed = rotateToA(patch, patchConvexity, a, b, c);
    
    a++;b++;c++;d++;


    // We'll construct the left rectangle

    // a nodes
        std::vector<int> anodes(a, 0);
        auto it = patch.begin();
        for (int i=0; i<a; i++){
            anodes[i] = Halfedge(m, *it).from();
            it++;
        }

    // b nodes
        int bsize = roundUpDivide(b, 2);
        if (b%2==0)
            bsize++;
        std::vector<int> bnodes(bsize, 0);
        it--;
        for (int i=0; i<bsize; i++){
            bnodes[i] = Halfedge(m, *it).from();
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

        createPointsBetween2Vx(cnodes, a-1, m, bvh);
        std::reverse(cnodes.begin(), cnodes.end());

    // Let's do d nodes now
        std::vector<int> dnodes(bsize, 0);
        dnodes[0]=Halfedge(m, *patch.begin()).from();
        it = patch.end();
        it--;
        for (int i=1; i<bsize; i++){
            dnodes[i] = Halfedge(m, *it).from();
            it--;
        }

    if (!wasReversed)
        meshingRectangle(anodes, bnodes, cnodes, dnodes, m, bvh);
    else
        meshingRectangle(dnodes, cnodes, bnodes, anodes, m, bvh);




    // Working on the right rectangle now. if the patch was too small (b<=2), we don't have a right triangle

    std::vector<int> anodes2(c, 0);
    if (b>2){
        // a nodes
        anodes = std::vector<int> (c,0);
        anodes2[0] = bnodes[bnodes.size()-1];
        it = patch.begin();
        std::advance(it, a+b+c-3+(b-1)/2);
        anodes2[c-1] = Halfedge(m, *it).from();

        createPointsBetween2Vx(anodes2, c-1, m, bvh);
        std::reverse(anodes2.begin(), anodes2.end());

        // b nodes
        std::vector bnodes2(roundUpDivide(b,2), 0);
        it = patch.begin();
        std::advance(it, a+b-2 - roundUpDivide(b,2)+1);
        for (int i=0; i<(int)bnodes2.size(); i++){
            bnodes2[i] = Halfedge(m, *it).from();
            it++;
        }

        // c nodes
        std::vector<int> cnodes2(c, 0);
        it--;
        for (int i=0; i<c; i++){
            cnodes2[i] = Halfedge(m, *it).from();
            it++;
        }
        std::reverse(cnodes2.begin(), cnodes2.end());

        // d nodes
        std::vector<int> dnodes2(bnodes2.size(), 0);
        it--;
        for (int i=0; i<(int)dnodes2.size(); i++){
            dnodes2[i] = Halfedge(m, *it).from();
            it++;
        }
        std::reverse(dnodes2.begin(), dnodes2.end());

        if (!wasReversed)
            meshingRectangle(anodes2, bnodes2, cnodes2, dnodes2, m, bvh);
        else
            meshingRectangle(dnodes2, cnodes2, bnodes2, anodes2, m, bvh);

    } else {

        // patch is too small to have a second rectangle, so here we only prepare the struture for the last triangle remesh
        it = patch.begin();
        std::advance(it, anodes.size()+bnodes.size()-2);
        for (int i=0; i<c; i++){
            anodes2[i] = Halfedge(m, *it).from();
            it++;
        }
        std::reverse(anodes2.begin(), anodes2.end());
    }


    // Remeshing the triangle part

    std::vector<int> btmPart(d-b+1, 0);
    it--;
    for (int i=0; i<(int)btmPart.size(); i++){
        btmPart[i] = Halfedge(m, *it).from();
        it++;
    }

    std::list<int> lst;
    std::reverse(anodes2.begin(), anodes2.end());
    lst.insert(lst.end(), anodes2.begin(), anodes2.end()-1);
    lst.insert(lst.end(), btmPart.begin(), btmPart.end()-1);
    lst.insert(lst.end(), cnodes.begin(), cnodes.end()-1);


    // Converting the vertices to halfedges for the nPatchRemesh function
    for (auto it = lst.begin(); it != lst.end(); ++it) {
        assert(Vertex(m, *it).halfedge().from() == Vertex(m, *it));
        *it = Vertex(m, *it).halfedge();
    }

    if (wasReversed){
        std::reverse(partSegments, partSegments + 6);
        std::reverse(std::next(lst.begin()), lst.end());
    }

    ajustPartSegments(partSegments, a-1, d-b, c-1);
    nPatchRemesh(partSegments, lst, m, 3, bvh);
}

inline bool remeshingPatch(std::list<int>& patch, std::list<int>& patchConvexity, int nEdge, Quads& m, FacetAttribute<int>& fa, int v, BVH bvh){
    assert(patchConvexity.front() >= 1);
    assert(nEdge == 3 || nEdge == 5 || nEdge == 4);

    // Segments contains the number of points between each edge of the patch
    // PartSegments are the segments but divided in 2 parts, according to the results of Bunin's equations
    int segments[] = {0,0,0,0,0};
    int partSegments[] = {0,0,0,0,0,0,0,0,0,0};
    segmentConstruction(patchConvexity, segments, nEdge);

    int a = 0;
    int b = 0;
    int c = 0;
    int d = 0;
    int solve4equationsCase = 0;

    if (nEdge == 4){
        solve4equationsCase = solve4equations(segments, partSegments, a, b, c, d);
        if (solve4equationsCase == 1){
            rectanglePatchRemesh(patch, segments, m, bvh);
            cleaningTopology(m, fa);
            std::cout << "solve " << nEdge << " (rectangle) equations success,    root: " << v << std::endl;
            return true;
        }

        if (solve4equationsCase == 2){
            quadrilateralPatchRemesh(partSegments, patch, patchConvexity, m, bvh, a, b, c, d);
            cleaningTopology(m, fa);
            std::cout << "solve " << nEdge << " (nonRectangle) equations success, root: " << v << std::endl;
            return true;
        }

    } else if (nEdge == 3){
        if (solve3equations(segments, partSegments)){
            nPatchRemesh(partSegments, patch, m, nEdge, bvh);
            cleaningTopology(m, fa);
            std::cout << "solve " << nEdge << " equations success,                root: " << v << std::endl;
            return true;
        }

    } else if (nEdge == 5){
        if (solve5equations(segments, partSegments)){
            nPatchRemesh(partSegments, patch, m, nEdge, bvh);
            cleaningTopology(m, fa);
            std::cout << "solve " << nEdge << " equations success,                root: " << v << std::endl;
            return true;
        }
    }

    return false;
}
