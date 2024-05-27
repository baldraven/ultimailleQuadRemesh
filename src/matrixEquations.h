#include <ultimaille/all.h>
// https://www.mcs.anl.gov/~fathom/meshkit-docs/html/Mesh_8cpp_source.html (Jaal)

inline int solve5equations(const int *segments, int *partSegments){
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

    // solved matrix computation
    std::vector<int> x(10);
    for (int i = 0; i < 10; i++) {
        double sum = 0.0;
        for (int j = 0; j < 10; j++)
            sum += M[i][j] * rhs[j];
        x[i] = (int) sum;
    }

    for (int i = 0; i < 10; i++){
        if (x[i] < 1) return 0;
    }

    if (x[0] != x[8]) return 0;
    if (x[0] + x[5] != rhs[5]) return 0;
    partSegments[0] = x[0];
    partSegments[1] = x[5];

    if (x[1] != x[9]) return 0;
    if (x[1] + x[6] != rhs[6]) return 0;
    partSegments[2] = x[1];
    partSegments[3] = x[6];

    if (x[2] != x[5]) return 0;
    if (x[2] + x[7] != rhs[7]) return 0;
    partSegments[4] = x[2];
    partSegments[5] = x[7];

    if (x[3] != x[6]) return 0;
    if (x[3] + x[8] != rhs[8]) return 0;
    partSegments[6] = x[3];
    partSegments[7] = x[8];

    if (x[4] != x[7]) return 0;
    if (x[4] + x[9] != rhs[9]) return 0;
    partSegments[8] = x[4];
    partSegments[9] = x[9];

    return 1;
}

inline int solve3equations(const int *segments, int *partsegments){
    double M[6][6], rhs[6];
    M[0][0] = -0.5;M[0][1] = -0.5;M[0][2] = 0.5;M[0][3] = 0.5;M[0][4] = 0.5;M[0][5] = -0.5;
    M[1][0] = 0.5;M[1][1] = -0.5;M[1][2] = -0.5;M[1][3] = -0.5;M[1][4] = 0.5;M[1][5] = 0.5;
    M[2][0] = -0.5;M[2][1] = 0.5;M[2][2] = -0.5;M[2][3] = 0.5;M[2][4] = -0.5;M[2][5] = 0.5;
    M[3][0] = 0.5;M[3][1] = 0.5;M[3][2] = -0.5;M[3][3] = 0.5;M[3][4] = -0.5;M[3][5] = 0.5;
    M[4][0] = -0.5;M[4][1] = 0.5;M[4][2] = 0.5;M[4][3] = 0.5;M[4][4] = 0.5;M[4][5] = -0.5;
    M[5][0] = 0.5;M[5][1] = -0.5;M[5][2] = 0.5;M[5][3] = -0.5;M[5][4] = 0.5;M[5][5] = 0.5;

    rhs[0] = 0.0;
    rhs[1] = 0.0;
    rhs[2] = 0.0;

    rhs[3] = segments[0];
    rhs[4] = segments[1];
    rhs[5] = segments[2];

    std::vector<int> x(6);
    for (int i = 0; i < 6; i++) {
        double sum = 0.0;
        for (int j = 0; j < 6; j++)
            sum += M[i][j] * rhs[j];
        x[i] = (int) sum;
    }

    for (int i = 0; i < 6; i++)
        if (x[i] < 1) return 0;

    if (x[0] + x[3] != rhs[3]) return 0;
    partsegments[0] = x[0];
    partsegments[1] = x[3];

    if (x[1] + x[4] != rhs[4]) return 0;
    partsegments[2] = x[1];
    partsegments[3] = x[4];

    if (x[2] + x[5] != rhs[5]) return 0;
    partsegments[4] = x[2];
    partsegments[5] = x[5];

    return 1;
}

inline int solve4equations(int* segments, int* partSegments, int &a, int &b, int &c, int &d){
    if (segments[0] == segments[2] && segments[1] == segments[3]){
        return 1;
    }
    a = fmax(segments[0], segments[2]);
    c = fmin(segments[0], segments[2]);
    b = fmin(segments[1], segments[3]);
    d = fmax(segments[1], segments[3]);

    int segmentsTri[] = {d-b,  c, a};
    if (solve3equations(segmentsTri, partSegments)){
        return 2;
    }

    std::swap(a, d);
    std::swap(b, c);

    // Sideway triangle insertion
    int segmentsTri2[] = {d-b,  c, a};
    if (solve3equations(segmentsTri2, partSegments)){
        std::cout << "alt_insertion" << std::endl;
        return 2;
    }

    return false;
}
