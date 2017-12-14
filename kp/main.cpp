#include <iostream>
#include <vector>
#include <assert.h>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <set>

using namespace std;

struct Point {
    int x, y, id;

    Point() {}

    Point(int x_, int y_) : x(x_), y(y_), id(-1) {}

    Point(int x_, int y_, int id_) : x(x_), y(y_), id(id_) {}

    //Note: id of + and - is undefined
    Point operator+(const Point &p) const {
        return Point(x + p.x, y + p.y);
    }

    Point operator-(const Point &p) const {
        return Point(x - p.x, y - p.y);
    }

    // += and -= preserves id
    Point &operator+=(const Point &p) {
        x += p.x;
        y += p.y;
        return *this;
    }

    Point &operator-=(const Point &p) {
        x -= p.y;
        y -= p.y;
        return *this;
    }

    long long operator*(const Point &b) const {
        return x * 1ll * b.x + y * 1ll * b.y;
    }

    // a % b > 0 if a is to the right of b
    long long operator%(const Point &p) const {
        return 1LL * x * p.y - 1LL * y * p.x;
    }

    bool operator==(const Point &p) const {
//        cout << x << " " << y << " " << id << endl;
//        cout << p.x << " " << p.y << " " << p.id << endl;
        int eqXY = make_pair(x, y) == make_pair(p.x, p.y);
        int eqId = id == p.id;
        assert(eqXY == eqId);

        return eqId;
    }

    bool operator!=(const Point &p) const {
        return !(*this == p);
    }

    double dist() const {
        return sqrt(pow((double) x, 2) + pow((double) y, 2));
    }

    double angle(const Point &p) const {
        auto dp = (double) ((*this) * p);
        dp /= dist();
        dp /= p.dist();
        return acos(dp);
    }
};

ostream &operator<<(ostream &stream, const Point &p) {
    stream << "(" << p.x << ", " << p.y << ", id=" << p.id << ")";
    return stream;
}

struct Edge {
    int fromId, toId;

    Edge() {}

    Edge(int fromId_, int toId_) : fromId(fromId_), toId(toId_) {}

    bool operator==(const Edge &e) const {
        return fromId == e.fromId && toId == e.toId;
    }
};

ostream &operator<<(ostream &stream, const Edge &e) {
    stream << "[edge from " << e.fromId << " to " << e.toId << "]";
    return stream;
}

struct Triangle {
    Edge edges[3];
};

bool is_right_turn(const Point &a, const Point &b, const Point &c) {
    return (b - a) % (c - a) > 0;
}

bool is_left_turn(const Point &a, const Point &b, const Point &c) {
    return (b - a) % (c - a) < 0;
}

typedef vector<Point> Polygon;


void print_polygon(Polygon &points) {
    printf("%d\n", (int) points.size());
    for (auto point : points) {
        printf("%d %d\n", point.x, point.y);
    }
}

struct PointCmp {
    bool operator()(const Point &a, const Point &b) {
        return a.x != b.x ? a.x < b.x : a.y < b.y;
    }
};

struct PointCmpByX {
    bool operator()(const Point &a, const Point &b) {
        return a.x < b.x;
    }
};

struct PointCmpByY {
    bool operator()(const Point &a, const Point &b) {
        return a.y < b.y;
    }
};

class Kirkpatrick {
private:
    int max_id;

public:
    Kirkpatrick(Polygon &p) {
//        // find max_id
//        max_id = max_element(p.begin(), p.end(), [](const Point &a, const Point &b) {
//            return a.id < b.id;
//        })->id;
//
//        // build triangulation
//        const auto triangulation = triangulate(p);
//        vector<Triangle> triangles = triangulation.first;
//        vector<Edge> edges = triangulation.second;
//
//        // add triangle's faces
//        for (auto &face : bound_polygon(p)) {
//            auto result = triangulate(face);
//            triangles.insert(triangles.end(), result.first.begin(), result.first.end());
//            edges.insert(edges.end(), result.second.begin(), result.second.end());
//        }
//
//        // TODO: build tree
    }

    bool contains(const Point &point) {
        // TODO: implement
        return false;
    }

    vector<Polygon> split_y(Polygon &P) {
        bound_polygon(P);

        unordered_map<int, Point> U; // U for universe
        for (auto &p: P) U[p.id] = p;

        struct PointHasher {
            size_t operator()(const Point &p) const {
                return hash<int>()(p.id);
            }
        };
        unordered_map<Point, int, PointHasher> pos;
        for (int i = 0; i < int(P.size()); i++) pos[P[i]] = i;
        struct EdgeHasher {
            size_t operator()(const Edge &e) const {
                return hash<int>()(e.fromId) ^ hash<int>()(e.toId);
            }
        };
        struct EdgeCompare {
            unordered_map<int, Point> &U;

            explicit EdgeCompare(unordered_map<int, Point> &U) : U(U) {}

            bool operator()(const Edge &el, const Edge &er) const {
                auto l = (U[el.toId].x + U[el.fromId].x) / 2., r = (U[er.toId].x + U[er.fromId].x) / 2.;
                return l < r;
            }
        };

        // vertex_type
        auto vertex_type = [P, &pos](const Point &p) {
            auto &r = P[(pos[p] + 1) % P.size()], &q = P[(pos[p] - 1 + P.size()) % P.size()];
            auto a = (q - p).angle(r - p);
            if (is_right_turn(p, q, r)) a = 2 * M_PI - a;
            if (r.y < p.y && q.y < p.y && a < M_PI) {
                return 0; // start
            } else if (r.y < p.y && q.y < p.y && a > M_PI) {
                return 1; // split
            } else if (r.y > p.y && q.y > p.y && a < M_PI) {
                return 2; // end
            } else if (r.y > p.y && q.y > p.y && a > M_PI) {
                return 3; // merge
            } else {
                return 4; // regular
            }
        };

        // T
        auto ecmp = [&U](const Edge &el, const Edge &er) {
            auto l = (U[el.toId].x + U[el.fromId].x) / 2., r = (U[er.toId].x + U[er.fromId].x) / 2.;
            return l < r;
        };
        set<Edge, decltype(ecmp)> T(ecmp);

        // Q
        auto pcmp = [](const Point &p, const Point &q) { return p.y < q.y || (p.y == q.y && p.x > q.x); };
        priority_queue<Point, vector<Point>, decltype(pcmp)> Q(pcmp, P);

        // D
        vector<Edge> D;

        // helper
        unordered_map<Edge, Point, EdgeHasher> helper;


        while (!Q.empty()) {
            auto vi = Q.top();
            Q.pop();

            auto type = vertex_type(vi);
            auto vi_1 = P[(pos[vi] - 1 + P.size()) % P.size()];  // v_{i - 1}
            auto vi1 = P[(pos[vi] + 1) % P.size()];  // v_{i + 1}
            auto ei = Edge(vi.id, vi1.id);  // e_i
            auto ei_1 = Edge(vi_1.id, vi.id);  // e_{i - 1}

//            cout << "vi_1 = " << vi_1 << endl;
//            cout << "vi = " << vi << endl;
//            cout << "vi1 = " << vi1 << endl;
//            cout << "type= " << type << endl;
//            cout << "ei_1= " << ei_1 << endl;
//            cout << "ei= " << ei << endl;
//            cout << "T size before= " << T.size() << endl;

            if (type == 0) { // start
                T.insert(ei);
                helper[ei] = vi;
            } else if (type == 1) { // split
                auto ej_it = --T.lower_bound(Edge(vi.id, vi.id));
                if (ej_it != T.end()) {
                    auto ej = *ej_it;
                    D.emplace_back(vi.id, helper[ej].id);
                    helper[ej] = vi;
                }
                T.insert(ei);
                helper[ei] = vi;
            } else if (type == 2) { // end
                if (helper.find(ei_1) != helper.end() && vertex_type(helper[ei_1]) == 3) {  // merge
                    D.emplace_back(vi.id, helper[ei_1].id);
                }
                T.erase(ei_1);
            } else if (type == 3) { // merge
                if (helper.find(ei_1) != helper.end() && vertex_type(helper[ei_1]) == 3) {  // merge
                    D.emplace_back(vi.id, helper[ei_1].id);
                }
                T.erase(ei_1);
                auto ej_it = --T.lower_bound(Edge(vi.id, vi.id));
                if (ej_it != T.end()) {
                    auto ej = *ej_it;
                    if (helper.find(ej) != helper.end() && vertex_type(helper[ej]) == 3) {  // merge
                        D.emplace_back(vi.id, helper[ej].id);
                    }
                    helper[ej] = vi;
                }
            } else { // regular
                bool cond = vi1.y < vi_1.y || (vi1.y == vi_1.y && vi1.x < vi_1.x);  // TODO: shitcases?
//                cout << "cond= " << cond << endl;
                if (cond) {
                    if (helper.find(ei_1) != helper.end() && vertex_type(helper[ei_1]) == 3) {  // merge
                        D.emplace_back(vi.id, helper[ei_1].id);
                    }
                    T.erase(ei_1);
                    T.insert(ei);
                    helper[ei] = vi;
                } else {
                    auto ej_it = --T.lower_bound(Edge(vi.id, vi.id));
                    if (ej_it != T.end()) {
                        auto ej = *ej_it;
                        if (helper.find(ej) != helper.end() && vertex_type(helper[ej]) == 3) {  // merge
                            D.emplace_back(vi.id, helper[ej].id);
                        }
                        helper[ej] = vi;
                    }
                }
            }

//            cout << "T size after= " << T.size() << endl;
//            cout << "------------------------" << endl;
        }

        // res
        vector<Polygon> res;
        int lastl, lastr;
        bool has_last_cut = false;
        for (auto &e : D) {
            auto l = pos[U[e.fromId]], r = pos[U[e.toId]];
            if (l > r) swap(l, r);
            if (has_last_cut) {
                vector<Point> v;
                v.insert(v.end(), P.begin() + l, P.begin() + lastl + 1);
                v.insert(v.end(), P.begin() + lastr, P.begin() + r + 1);
                res.push_back(v);
            } else {
                res.emplace_back(vector<Point>(P.begin() + l, P.begin() + r + 1));
            }
            lastl = l;
            lastr = r;
            has_last_cut = true;
        }
        if (has_last_cut) {
            vector<Point> v;
            v.insert(v.end(), P.begin(), P.begin() + lastl);
            v.insert(v.end(), P.begin() + lastr, P.end());
            v.push_back(P[0]);
            res.push_back(v);
        } else {
            res.push_back(P);
        }
        return res;
    }

    pair<vector<Triangle>, vector<Edge>> triangulate_y(Polygon points) {
        // TODO: paste Seva's implementation
        return {{},
                {}};
    }

    pair<vector<Triangle>, vector<Edge>> triangulate(Polygon points) {
        vector<Triangle> triangles;
        vector<Edge> edges;

        for (auto &face : split_y(points)) {
            const auto result = triangulate_y(face);
            triangles.insert(triangles.end(), result.first.begin(), result.first.end());
            edges.insert(edges.end(), result.second.begin(), result.second.end());
        }

        return {triangles, edges};
    }

    //orders p in counter-clockwise order. works for non-convex polygons
    void reorder_ccw(Polygon &p) {
        assert(p.size() >= 3);

        long long square = 0;
        for (int i = 1; i < (int) p.size() - 1; i++) {
            square += (p[i] - p[0]) % (p[i + 1] - p[0]);
        }

        if (square < 0) {
            reverse(p.begin(), p.end());
        }
    }

    // find additional faces
    vector<Polygon> bound_polygon(Polygon &points) {
        reorder_ccw(points);
        rotate(points.begin(), min_element(points.begin(), points.end(), PointCmp()), points.end());

        // build convex hull
        const auto hull = convex_hull(points);

        // find left bottom corner
        const int min_x = min_element(hull.begin(), hull.end(), PointCmpByX())->x;
        const int min_y = min_element(hull.begin(), hull.end(), PointCmpByY())->y;

        // find the farthest point
        auto farthest = *max_element(hull.begin(), hull.end(), [](const Point &a, const Point &b) {
            return a.x + a.y < b.x + b.y;
        });

        // find other coordinates
        const int max_x = farthest.x + (farthest.y - min_y);
        const int max_y = farthest.y + (farthest.x - min_x);

        const Polygon triangle = {Point(min_x, min_y, ++max_id),
                                  Point(max_x, min_y, ++max_id),
                                  Point(min_x, max_y, ++max_id)};

        int state = 0;
        int iter = 0;

        vector<Polygon> faces;

        for (int i = 0; i < hull.size(); ++i) {
            int j = (i == hull.size() - 1 ? 0 : i + 1);

            Polygon face;

            while (points[iter] != hull[j]) {
                face.push_back(points[iter]);
                iter = (iter == points.size() - 1 ? 0 : iter + 1);
            }

            face.push_back(points[iter]);

            bool is_vertical = (hull[j].x == hull[i].x && hull[j].y < hull[i].y);
            bool is_horizontal = (hull[j].x > hull[i].x && hull[j].y == hull[i].y);
            bool is_45_degree = (hull[j].x < hull[i].x && hull[i].x - hull[j].x == hull[j].y - hull[i].y);
            bool is_any = is_vertical || is_horizontal || is_45_degree;

            if (is_any) {
                state = (is_vertical ? 0 : (is_horizontal ? 1 : 2));
            } else {
                face.push_back(triangle[state]);
            }

            if (faces.size() >= 3) {
                faces.push_back(face);
            }
        }

        return faces;
    }

    Polygon convex_hull(Polygon points) {
        sort(points.begin(), points.end(), PointCmp());

        if (points.size() <= 3) {
            return points;
        }

        vector<Point> upper = {points.front()};
        vector<Point> lower = {points.front()};

        for (int i = 1; i < points.size(); ++i) {
            if (i == points.size() - 1 || is_right_turn(points.front(), points[i], points.back())) {
                while (upper.size() > 1 && is_left_turn(upper[upper.size() - 2], upper.back(), points[i])) {
                    upper.pop_back();
                }

                upper.push_back(points[i]);
            }

            if (i == points.size() - 1 || is_left_turn(points.front(), points[i], points.back())) {
                while (lower.size() > 1 && is_right_turn(lower[lower.size() - 2], lower.back(), points[i])) {
                    lower.pop_back();
                }

                lower.push_back(points[i]);
            }
        }

        for (int i = (int) lower.size() - 2; i > 0; --i) {
            upper.push_back(lower[i]);
        }

        return upper;
    }
};

int main() {
    // test cases
    Polygon P1 = {{2, 0, 0},
                  {3, 2, 1},
                  {6, 0, 2},
                  {2, 5, 3},
                  {2, 3, 4},
                  {0, 3, 5}};
    Polygon P2 = {{0, 1, 0},
                  {1, 0, 1},
                  {1, 1, 2}};
    auto P = P1;

    auto k = Kirkpatrick(P);
    auto res = k.split_y(P);

    for (int i = 0; i < int(res.size()); i++) {
        cout << i + 1 << "th polygon:" << endl;
        print_polygon(res[i]);
        if (i != int(res.size()) - 1) {
            cout << "------------" << endl;
        }
    }
}