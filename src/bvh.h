#pragma once

#include <ultimaille/all.h>
#include <numeric>
#include <algorithm>


// From Yoann Coudert-Osmont

using namespace UM;

struct BVH {
    using Box = std::array<double, 6>;

private:
    const UM::Triangles &m;
    std::vector<std::tuple<int, int, Box>> nodes;

    inline vec3 geom(const int h) { return m.points[m.vert(h/3, (h+1)%3)] - m.points[m.facets[h]]; }

public:
    inline Box tri_box(int f) {
        return {
            std::min({m.points[m.vert(f, 0)][0], m.points[m.vert(f, 1)][0], m.points[m.vert(f, 2)][0]}),
            std::max({m.points[m.vert(f, 0)][0], m.points[m.vert(f, 1)][0], m.points[m.vert(f, 2)][0]}),
            std::min({m.points[m.vert(f, 0)][1], m.points[m.vert(f, 1)][1], m.points[m.vert(f, 2)][1]}),
            std::max({m.points[m.vert(f, 0)][1], m.points[m.vert(f, 1)][1], m.points[m.vert(f, 2)][1]}),
            std::min({m.points[m.vert(f, 0)][2], m.points[m.vert(f, 1)][2], m.points[m.vert(f, 2)][2]}),
            std::max({m.points[m.vert(f, 0)][2], m.points[m.vert(f, 1)][2], m.points[m.vert(f, 2)][2]})
        };
    }

    inline static void surround(Box &a, const Box &b) {
        for(int i = 0; i < 6; i += 2) a[i] = std::min(a[i], b[i]);
        for(int i = 1; i < 6; i += 2) a[i] = std::max(a[i], b[i]);
    }

    inline static std::pair<double, double> volume_box(const Box &box) {
        return {(box[1]-box[0]) * (box[3]-box[2]) * (box[5]-box[4]), 
                std::max({ (box[1]-box[0]) * (box[3]-box[2]), (box[1]-box[0]) * (box[5]-box[4]), (box[5]-box[4]) * (box[3]-box[2]) })};
    }

    inline double dist_segment(double a, double b, double x) { return x < a ? a-x : (x > b ? x-b : 0.); }
    inline double dist2_box(const Box &box, const vec3 &p) {
        return vec3(dist_segment(box[0], box[1], p.x), dist_segment(box[2], box[3], p.y), dist_segment(box[4], box[5], p.z)).norm2();
    }

    vec3 proj_facet(const vec3 &p, int f) {
        const vec3 n = UM::cross(geom(3*f+0), geom(3*f+1)).normalized();
        vec3 q = p - ((p-m.points[m.vert(f,0)])*n) * n;
        for(int i = 0; i < 3; ++i) {
            const vec3 ne = UM::cross(geom(3*f+i), n).normalized();
            const double d = ne * (q - m.points[m.vert(f,i)]);
            if(d > 0.) q -= d*ne;
        }
        return q;
    }

    void init(std::vector<int>::iterator start, size_t n) {
        const size_t ind = nodes.size();
        nodes.emplace_back();
        int axis;
        const auto maxAxis = [&](int f)->double {
            return std::max({m.points[m.vert(f, 0)][axis], m.points[m.vert(f, 1)][axis], m.points[m.vert(f, 2)][axis]});
        };
        const auto comp = [&](int a, int b) {
            return maxAxis(a) < maxAxis(b);
        };

        int bestAxis = 0;
        int bestSep = 1;
        std::pair<double, double> bestScore = { std::numeric_limits<double>::max(), std::numeric_limits<double>::max() };
        std::vector<std::pair<double, double>> volumes(n-1);
        Box &box = std::get<2>(nodes[ind]);
        for(axis = 0; axis < 3; ++axis) {
            sort(start, start+n, comp);
            int i = 0;
            box = tri_box(*start);
            while(true) {
                volumes[i] = volume_box(box);
                if(++i < (int)volumes.size()) surround(box, tri_box(*(start+i)));
                else break;
            }
            box = tri_box(*(start+i));
            while(true) {
                std::pair<double, double> score = volume_box(box);
                score.first = score.first * (n - i) + volumes[i-1].first * i;
                score.second = score.second * (n - i) + volumes[i-1].second * i;
                if(score < bestScore) {
                    bestScore = score;
                    bestAxis = axis;
                    bestSep = i;
                }
                if(--i > 0) surround(box, tri_box(*(start+i)));
                else break;
            }
        }

        if(bestAxis < 2) {
            axis = bestAxis;
            std::sort(start, start+n, comp);
        }

        std::get<0>(nodes[ind]) = nodes.size();
        if(bestSep == 1) nodes.emplace_back(*start, -1, tri_box(*start));
        else init(start, bestSep);
        std::get<1>(nodes[ind]) = nodes.size();
        if(bestSep+1 == (int)n) nodes.emplace_back(*(start + bestSep), -1, tri_box(*(start + bestSep)));
        else init(start + bestSep, n-bestSep);
        std::get<2>(nodes[ind]) = std::get<2>(nodes[std::get<0>(nodes[ind])]);
        surround(std::get<2>(nodes[ind]), std::get<2>(nodes[std::get<1>(nodes[ind])]));
    }

    BVH(const UM::Triangles &m): m(m), nodes() {
        if(m.nfacets() < 2) return;
        std::vector<int> tri_inds(m.nfacets());
        std::iota(tri_inds.begin(), tri_inds.end(), 0);
        init(tri_inds.begin(), m.nfacets());
    }

    vec3 project(const vec3 &p) {
        if(m.nfacets() == 0) return p;
        else if(m.nfacets() == 1) return proj_facet(p, 0);
        using QEl = std::pair<double, int>;
        std::priority_queue<QEl, std::vector<QEl>, std::greater<QEl>> Q;
        Q.emplace(0., 0);
        double best_dist2 = std::numeric_limits<double>::max();
        vec3 q;
        while(!Q.empty()) {
            if(Q.top().first >= best_dist2) break;
            const int i = Q.top().second; Q.pop();
            for(int j : {std::get<0>(nodes[i]), std::get<1>(nodes[i])}) {
                if(std::get<1>(nodes[j]) == -1) {
                    const vec3 q2 = proj_facet(p, std::get<0>(nodes[j]));
                    const double d2 = (p-q2).norm2();
                    if(d2 < best_dist2) {
                        best_dist2 = d2;
                        q = q2;
                    }
                } else Q.emplace(dist2_box(std::get<2>(nodes[j]), p), j); 
            }
        }
        return q;
    }
};