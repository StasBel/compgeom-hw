#include <functional>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cassert>
#include <queue>
#include <cmath>
#include <unordered_map>
#include <set>
#include <unordered_set>

#pragma comment(linker, "/STACK: 2000000")

using namespace std;

struct Point {
    int x, y, id;

    Point() {}

    Point(int x_, int y_) : x(x_), y(y_), id(-1) {}

    Point(int x_, int y_, int id_) : x(x_), y(y_), id(id_) {}

    //Note: _id of + and - is undefined
    Point operator+(const Point &p) const {
        return Point(x + p.x, y + p.y);
    }

    Point operator-(const Point &p) const {
        return Point(x - p.x, y - p.y);
    }

    // += and -= preserves _id
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

    double fatan() const {
        auto a = atan2((double) y, (double) (x == 0 && y == 0) ? 1 : x);
//        auto a = atan2(y, x);
        if (a < 0) a += 2 * M_PI;
        return a;
    }

    double atan() const {
        return atan2((double) y, (double) (x == 0 && y == 0) ? 1 : x);
    }
};

ostream &operator<<(ostream &stream, const Point &p) {
    stream << "(" << p.x << ", " << p.y << ", _id=" << p.id << ")";
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
        printf("%d %d _id=%d\n", point.x, point.y, point.id);
    }
}

inline Point getPoint(const Polygon &p, int i) {
    i = i < 0 ? i + (int) p.size() : i;
    i = i >= (int) p.size() ? i - (int) p.size() : i;
    return p[i];
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

inline int sign(long long x) {
    return x == 0 ? 0 : (x < 0 ? -1 : 1);
}

class Kirkpatrick {
private:
    class Node {
    public:
        Node(Triangle triangle, size_t id)
                : _triangle(triangle), _children(), _id(id) {}

        Node(Triangle triangle, vector<size_t> children, size_t id)
                : _triangle(triangle), _children(children), _id(id) {}

        Triangle getTriangle() const {
            return _triangle;
        }

        vector<size_t> getChildren() const {
            return _children;
        }

        void addChild(size_t nodeId) {
            _children.push_back(nodeId);
        }

        size_t getId() const {
            return _id;
        }

    private:
        Triangle _triangle;
        vector<size_t> _children;
        size_t _id;
    };

    class Vertex {
    public:

        explicit Vertex(const Point &point)
                : _deg(0), _point(point), _incidentEdges(), _incidentTriangles() {}

        Vertex() {}

        int _deg;
        Point _point;
        vector<size_t> _incidentEdges;
        vector<size_t> _incidentTriangles;

        Point &getPoint() {
            return _point;
        }

        int getId() const {
            return _point.id;
        }

        void decDegree() {
            --_deg;
        }

        vector<size_t> &incidentEdges() {
            return _incidentEdges;
        }

        vector<size_t> &incidentTriangles() {
            return _incidentTriangles;
        }

        void addInEdge(size_t edgeIndex) {
            _incidentEdges.push_back(edgeIndex);
        }

        void addOutEdge(size_t edgeIndex) {
            _incidentEdges.push_back(edgeIndex);
            ++_deg;
        }

        void addTriangle(size_t triangleIndex) {
            if (!_incidentTriangles.empty() && _incidentTriangles.back() == triangleIndex) {
                return;
            }
            _incidentTriangles.push_back(triangleIndex);
        }

        int degree() const {
            return _deg;
        }
    };

private:
    int max_id;
    Polygon polygon;
    vector<Node> dgraph;
    static const int K = 6;
    size_t outsideTrianglesIdBegin;
    size_t outsideTrianglesIdEnd;


    bool isBoundTrianglePoint(int pointId) {
        return max_id - 3 < pointId && pointId <= max_id;
    }

    bool isIntersected(const Triangle &first, const Triangle &second) {
        bool result = false;
        for_each(first.edges, first.edges + 3, [this, &result, &second](const Edge &firstEdge) {
            bool intersect = triangleContains(second, polygon[firstEdge.fromId]) ||
                             triangleContains(second, polygon[firstEdge.toId]);
            if (intersect) {
                result = true;
            }
        });
        return result;
    }

    bool triangleContains(const Triangle &triangle, const Point &point) {
        bool result = true;
        for_each(triangle.edges, triangle.edges + 3, [this, &result, &point](const Edge &edge) {
            result &= isLeftRotation(polygon[edge.fromId], polygon[edge.toId], point);
        });
        return result;
    }

    bool isLeftRotation(const Point &p, const Point &q, const Point &r) {
        int64_t temp = 1ll * (q.x - p.x) * (r.y - q.y) - 1ll * (r.x - q.x) * (q.y - p.y);
        return temp >= 0;
    }

    Edge getNonIncidentEdge(const Triangle &triangle, int vertexId) {
        for (auto edge : triangle.edges) {
            if (edge.fromId != vertexId && edge.toId != vertexId) {
                return edge;
            }
        }
        assert(false);
    }

public:
    Kirkpatrick(Polygon &p)
            : polygon(p) {
        // find max_id
        max_id = max_element(p.begin(), p.end(), [](const Point &a, const Point &b) {
            return a.id < b.id;
        })->id;

        // build triangulation
        const auto triangulation = triangulate(p);
        vector<Triangle> triangles = triangulation.first;
        vector<Edge> edges = triangulation.second;

        // using in contains, to check whether current _triangle inside polygon or not
        outsideTrianglesIdBegin = triangles.size();

        // add _triangle's faces
        for (auto &face : bound_polygon(p)) {
            auto result = triangulate(face);
            triangles.insert(triangles.end(), result.first.begin(), result.first.end());
            edges.insert(edges.end(), result.second.begin(), result.second.end());
        }

        // also using in contains
        outsideTrianglesIdEnd = triangles.size();
        buildDGraph(triangles, edges);
    }

    void buildDGraph(vector<Triangle> triangles, vector<Edge> edges) {
        for_each(triangles.begin(), triangles.end(), [this](const Triangle &triangle) {
            size_t id = dgraph.size();
            dgraph.push_back(Node(triangle, id));
        });

        unordered_map<int, Vertex> idToVertex;

        for_each(polygon.begin(), polygon.end(), [&idToVertex](const Point &point) {
            idToVertex[point.id] = Vertex(point);
        });

        for (size_t i = 0; i < edges.size(); ++i) {
            idToVertex[edges[i].fromId].addOutEdge(i);
            idToVertex[edges[i].toId].addInEdge(i);
        }

        for (size_t i = 0; i < triangles.size(); ++i) {
            for_each(triangles[i].edges, triangles[i].edges + 3, [&idToVertex, &i](const Edge &edge) {
                idToVertex[edge.fromId].addTriangle(i);
            });
        }

        vector<bool> deletedTriangles(triangles.size(), false);
        vector<bool> deletedEdges(edges.size(), false);
        while (idToVertex.size() > 3) {
            unordered_set<int> incidentVertex{};
            vector<int> deletedVertex{};
            for (std::pair<const int, Vertex> &vertexIt: idToVertex) {
                int vertexId = vertexIt.first;
                Vertex vertex = vertexIt.second;
                if (incidentVertex.find(vertexId) == incidentVertex.end() &&
                    !isBoundTrianglePoint(vertexId) && vertex.degree() < K) {

                    for_each(vertex.incidentEdges().begin(), vertex.incidentEdges().end(),
                             [&edges, &idToVertex, &deletedEdges, &incidentVertex](size_t edgeIndex) {
                                 if (!deletedEdges[edgeIndex]) {
                                     incidentVertex.insert(edges[edgeIndex].fromId);
                                     incidentVertex.insert(edges[edgeIndex].toId);
                                     idToVertex[edges[edgeIndex].fromId].decDegree();
                                     deletedEdges[edgeIndex] = true;
                                 }
                             }
                    );
//                    std::cout << "delete vertex " << vertexId << std::endl;
                    vector<Edge> polygonSide{};
                    Polygon polygon4Trian{};
                    vector<size_t> deletedTrianglesId{};

                    for_each(vertex.incidentTriangles().begin(),
                             vertex.incidentTriangles().end(),
                             [this, &idToVertex, &polygonSide,
                                     &deletedTriangles, &deletedTrianglesId, &vertexId](size_t triangleId) {
                                 if (!deletedTriangles[triangleId]) {
                                     Edge nonIncidentEdge = getNonIncidentEdge(dgraph[triangleId].getTriangle(),
                                                                               vertexId);
                                     polygonSide.push_back(nonIncidentEdge);
                                     idToVertex[nonIncidentEdge.fromId].decDegree();
                                     deletedTrianglesId.push_back(triangleId);
                                     deletedTriangles[triangleId] = true;
                                 }
                             }
                    );

                    Point point = polygon[vertex.getId()];
                    std::sort(polygonSide.begin(), polygonSide.end(),
                              [this, &idToVertex, &point](const Edge &first, const Edge &second) {
                                  Point a = idToVertex[first.fromId].getPoint();
                                  Point b = idToVertex[second.fromId].getPoint();

                                  if (a.x - point.x >= 0 && b.x - point.x < 0) {
                                      return false;
                                  }
                                  if (a.x - point.x < 0 && b.x - point.x >= 0) {
                                      return true;
                                  }
                                  if (a.x - point.x == 0 && b.x - point.x == 0) {
                                      return a.y < b.y;
                                  }

                                  int64_t det = 1ll * (a.x - point.x) * (b.y - point.y) -
                                                1ll * (b.x - point.x) * (a.y - point.y);
                                  if (det > 0) {
                                      return true;
                                  }
                                  if (det < 0) {
                                      return false;
                                  }
                                  int64_t dist1 =
                                          1ll * (a.x - point.x) * (a.x - point.x) + (a.y - point.y) * (a.y - point.y);
                                  int64_t dist2 =
                                          1ll * (b.x - point.x) * (b.x - point.x) + (b.y - point.y) * (b.y - point.y);

                                  return dist1 < dist2;
                              });
                    for_each(polygonSide.begin(), polygonSide.end(), [this, &polygon4Trian](const Edge &edge) {
                        polygon4Trian.push_back(polygon[edge.fromId]);
                    });
//                    std::cout << "point for polygon\n";
//                    for_each(polygon4Trian.begin(), polygon4Trian.end(), [] (const Point& point) {
//                        std::cout << point.x << "  " << point.y << std::endl;
//                    });

                    auto triangulation = triangulate(polygon4Trian);
                    vector<Triangle> insideTriangles = triangulation.first;
                    vector<Edge> insideEdges = triangulation.second;

//                    std::cout << "getting triangulation\n";
//                    for_each(insideTriangles.begin(), insideTriangles.end(), [this] (const Triangle& triangle) {
//                        for_each(triangle.edges, triangle.edges + 3, [this] (const Edge& edge) {
//                            std::cout << polygon[edge.fromId].x << " " << polygon[edge.fromId].y << " ";
//                        });
//                        std::cout << "\n";
//                    });

                    for_each(insideTriangles.begin(), insideTriangles.end(),
                             [this, &triangles, &idToVertex,
                                     &deletedTrianglesId, &deletedTriangles](Triangle &triangle) {
                                 size_t triangleId = dgraph.size();
                                 deletedTriangles.push_back(false);
                                 dgraph.push_back(Node(triangle, triangleId));

                                 for_each(triangle.edges, triangle.edges + 3,
                                          [&triangleId, &idToVertex](Edge &edge) {
                                              idToVertex[edge.fromId].addTriangle(triangleId);
                                          }
                                 );
                                 for_each(deletedTrianglesId.begin(), deletedTrianglesId.end(),
                                          [this, &triangleId, &triangles](size_t deletedTriangleId) {
                                              if (isIntersected(dgraph[triangleId].getTriangle(),
                                                                dgraph[deletedTriangleId].getTriangle())) {

                                                  dgraph[triangleId].addChild(deletedTriangleId);
                                              }
                                          }
                                 );

                             }
                    );
                    for_each(insideEdges.begin(), insideEdges.end(),
                             [&idToVertex, &edges, &deletedEdges](Edge &edge) {
                                 size_t edgeIndex = edges.size();
                                 edges.push_back(edge);
                                 deletedEdges.push_back(false);
                                 idToVertex[edge.fromId].addOutEdge(edgeIndex);
                                 idToVertex[edge.toId].addInEdge(edgeIndex);
                             }
                    );
                    deletedVertex.push_back(vertexId);
                }
            }
            for_each(deletedVertex.begin(), deletedVertex.end(), [&idToVertex](int id) {
                idToVertex.erase(id);
            });
        }
    }

    bool contains(const Point &point) {
        Node root = dgraph.back();
        while (!root.getChildren().empty()) {
            if (!triangleContains(root.getTriangle(), point)) {
                return false;
            }
            for (size_t childNodeId: root.getChildren()) {
                Node childNode = dgraph[childNodeId];
                if (triangleContains(childNode.getTriangle(), point)) {
                    root = childNode;
                    break;
                }
            }
        }
        return !(root.getId() >= outsideTrianglesIdBegin && root.getId() <= outsideTrianglesIdEnd);
    }

    vector<Polygon> split_y(Polygon &P) {
        reorder_ccw(P);

        unordered_map<int, Point> U;  // U is for universe
        for (auto &p: P) U[p.id] = p;

        unordered_map<int, int> pos;
        for (int i = 0; i < int(P.size()); i++) pos[P[i].id] = i;
        struct EdgeHasher {
            size_t operator()(const Edge &e) const {
                return hash<int>()(e.fromId) ^ hash<int>()(e.toId);
            }
        };

        // T
        auto pcmp = [](const Point &p, const Point &q) { return p.y < q.y || (p.y == q.y && p.x > q.x); };
        auto ecmp = [&](const Edge &el, const Edge &er) {
            auto a = U[el.fromId], b = U[el.toId];
            if (pcmp(b, a)) swap(a, b);
            auto c = U[er.fromId], d = U[er.toId];
            if (pcmp(d, c)) swap(c, d);
            assert(a != b || c != d);
            if (a == b) {
                return is_right_turn(c, d, a);
            } else if (c == d) {
                return is_left_turn(a, b, c);
            } else {
                auto f = is_right_turn(c, d, a) && is_right_turn(c, d, b);
                auto s = is_left_turn(a, b, c) && is_left_turn(a, b, d);
                return s || f;
            }
        };
        set<Edge, decltype(ecmp)> T(ecmp);

        // Q
        priority_queue<Point, vector<Point>, decltype(pcmp)> Q(pcmp, P);

        // vertex_type
        auto vertex_type = [&](const Point &p) {
            auto &r = P[(pos[p.id] + 1) % P.size()], &q = P[(pos[p.id] - 1 + P.size()) % P.size()];
            auto a = (q - p).angle(r - p);
            if (is_right_turn(p, q, r)) a = 2 * M_PI - a;
            if (pcmp(q, p) && pcmp(r, p) && a < M_PI) {
                return 0; // start
            } else if (pcmp(q, p) && pcmp(r, p) && a > M_PI) {
                return 1; // split
            } else if (pcmp(p, q) && pcmp(p, r) && a < M_PI) {
                return 2; // end
            } else if (pcmp(p, q) && pcmp(p, r) && a > M_PI) {
                return 3; // merge
            } else {
                return 4; // regular
            }
        };

        // D
        vector<Edge> D;

        // helper
        unordered_map<Edge, Point, EdgeHasher> helper;

        while (!Q.empty()) {
            auto vi = Q.top();
            Q.pop();

            auto type = vertex_type(vi);
            auto vi_1 = P[(pos[vi.id] - 1 + P.size()) % P.size()];  // v_{i - 1}
            auto vi1 = P[(pos[vi.id] + 1) % P.size()];  // v_{i + 1}
            auto ei = Edge(vi.id, vi1.id);  // e_i
            auto ei_1 = Edge(vi_1.id, vi.id);  // e_{i - 1}

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
        }

//        // D Output
//        cout << "D.size() is " << D.size() << endl;
//        cout << D.size() << endl;
//        for (auto &e: D) {
//            cout << P[pos[e.fromId]].x << " " << P[pos[e.fromId]].y << " ";
//            cout << P[pos[e.toId]].x << " " << P[pos[e.toId]].y << endl;
//        }

        // construct graph
        int cur_base = -1;
        auto cmp = [&](int ll, int rr) {
            auto &p = U[cur_base], &q = U[ll], &r = U[rr];
            auto p1 = (q - p), p2 = (r - p);
            return p1.fatan() > p2.fatan();
        };
        vector<set<int, decltype(cmp)>> G(P.size(), set<int, decltype(cmp)>(cmp));
        for (int i = 0; i + 1 < int(P.size()); i++) {
            auto l = P[i].id, r = P[i + 1].id;
            cur_base = l;
            G[pos[l]].insert(r);
            cur_base = r;
            G[pos[r]].insert(l);
        }
        auto l = P[0].id, r = P[int(P.size()) - 1].id;
        cur_base = l;
        G[pos[l]].insert(r);
        cur_base = r;
        G[pos[r]].insert(l);
        for (auto &e : D) {
            auto l = U[e.fromId].id, r = U[e.toId].id;
            cur_base = l;
            G[pos[l]].insert(r);
            cur_base = r;
            G[pos[r]].insert(l);
        }

        // run dfs to collect faces
        vector<vector<int>> resi;
        int start = -1;
        vector<int> cur_pts;
        function<void(int, int)> dfs = [&](int x, int p) {
            cur_base = x;
            if (x == start) {
                resi.push_back(cur_pts);
                cur_pts.clear();
                start = x;
                cur_pts.push_back(x);
            } else {
                if (start == -1) start = x;
                cur_pts.push_back(x);
            }
            int px = pos[x];
            if (!G[px].empty()) {
                auto it = G[px].upper_bound(p);
                if (it == G[px].end()) it = G[px].begin();
                int n = *it;
                G[px].erase(n);
                dfs(n, x);
            }
        };
        for (auto &p : P) {
            cur_base = p.id;
            while (!G[pos[p.id]].empty()) {
                start = -1;
                cur_pts.clear();
                dfs(p.id, p.id);
                cur_base = p.id;
            }
        }

        // delete outer face and map id to points
        vector<Polygon> res;
        bool meet_outer = false;
        for (auto &pi: resi) {
            if (!meet_outer && pi.size() == P.size()) {
                meet_outer = true;
                continue;
            }
            Polygon nP;
            nP.reserve(pi.size());
            for (auto i: pi) {
                nP.push_back(U[i]);
            }
            assert(nP.size() >= 3);
            reorder_ccw(nP);
            res.push_back(nP);
        }
        assert(meet_outer);
        return res;
    }

    pair<vector<Triangle>, vector<Edge>> triangulate_y(Polygon points) {
        assert(points.size() >= 3);
        reorder_ccw(points);

        unordered_map<int, Point> nextPoint;
        for (int i = 0; i < (int) points.size(); i++) {
            nextPoint[points[i].id] = getPoint(points, i + 1);
        }

        order_polygon(points);

        vector<Triangle> triangles;
        vector<Edge> edges;
        vector<Point> curStack;
        for (const Point &pnt: points) {
            if (curStack.size() < 2) {
                curStack.push_back(pnt);
                continue;
            }
            Point last = curStack.back();
            curStack.pop_back();

            Point top = curStack[0];
//        cerr << pnt._id << " " << top._id << " " << nextPoint[pnt._id]._id << " " << nextPoint[top._id]._id << endl;
            if (nextPoint[pnt.id] == top || nextPoint[top.id] == pnt) {
                for (int i = 0; i < (int) curStack.size() - 1; i++) {
                    add_triangle(triangles, edges, pnt, curStack[i], curStack[i + 1]);
                }
                add_triangle(triangles, edges, pnt, curStack.back(), last);
                curStack.clear();

                curStack.push_back(last);
            } else {
                int expectedSign = last == nextPoint[pnt.id] ? 1 : -1;
                Point cur = last;
                while (curStack.size() > 0) {
                    if (sign((cur - pnt) % (curStack.back() - pnt)) != expectedSign) {
                        break;
                    }
//                    cerr << cur.id << " " << curStack.back().id << " " << pnt.id << " " << expectedSign << endl;
                    add_triangle(triangles, edges, cur, curStack.back(), pnt);
                    cur = curStack.back();
                    curStack.pop_back();
                }

                curStack.push_back(cur);
            }
            curStack.push_back(pnt);
            //       cerr << curStack.size() << endl;
        }

        //   cerr << curStack.size() << endl;
        assert(curStack.size() == 2);
        return {triangles, edges};
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

        // find left bottom corner
        const int min_x = min_element(points.begin(), points.end(), PointCmpByX())->x - 1;
        const int min_y = min_element(points.begin(), points.end(), PointCmpByY())->y - 1;

        // find the farthest point
        auto farthest = *max_element(points.begin(), points.end(), [](const Point &a, const Point &b) {
            return a.x + a.y < b.x + b.y;
        });

        // find other coordinates
        const int max_x = farthest.x + (farthest.y - min_y) + 1;
        const int max_y = farthest.y + (farthest.x - min_x) + 1;

        const Polygon triangle = {Point(min_x, min_y, ++max_id),
                                  Point(max_x, min_y, ++max_id),
                                  Point(min_x, max_y, ++max_id)};

        polygon.insert(polygon.end(), triangle.begin(), triangle.end());

        int iter = 0;

        Polygon p1;
        p1.push_back(triangle[0]);
        while (iter == 0 || (iter < points.size() && points[iter].id != farthest.id)) {
            p1.push_back(points[iter++]);
        }

        if (iter < points.size()) {
            p1.push_back(points[iter]);
        }

        p1.push_back(triangle[1]);

        Polygon p2;
        p2.push_back(triangle[1]);
        p2.push_back(farthest);
        p2.push_back(triangle[2]);

        Polygon p3;
        p3.push_back(triangle[2]);
        while (iter < points.size()) {
            p3.push_back(points[iter++]);
        }
        p3.push_back(points[0]);
        p3.push_back(triangle[0]);

        vector<Polygon> faces = {p1, p2, p3};
        return faces;
    }

    void add_triangle(vector<Triangle> &triangles, vector<Edge> &edges, Point p1, Point p2, Point p3) {
        long long prod = (p2 - p1) % (p3 - p2);

        assert(prod != 0);

        Triangle t;
        Point p[4] = {p1, p2, p3};

        if (prod < 0) {
            reverse(p, p + 3);
        }
        p[3] = p[0];

        for (int i = 0; i < 3; i++) {
            Edge e = Edge(p[i].id, p[i + 1].id);
            t.edges[i] = e;
            edges.push_back(e);
        }

        triangles.push_back(t);
    }

    void order_polygon(Polygon &points) {
        const auto &cmpByY = [](const Point &p1, const Point &p2) { return p1.y > p2.y; };
        // min, max -- not an error
        int maxPointPos = min_element(points.begin(), points.end(), cmpByY) - points.begin();
        int minPointPos = max_element(points.begin(), points.end(), cmpByY) - points.begin();

        vector<Point> path[2];
        for (int i = minPointPos; i != maxPointPos; i = (i + 1) % points.size()) {
            path[0].push_back(points[i]);
        }

        for (int i = maxPointPos; i != minPointPos; i = (i + 1) % points.size()) {
            path[1].push_back(points[i]);
        }

        reverse(path[0].begin(), path[0].end());

        merge(path[0].begin(), path[0].end(), path[1].begin(), path[1].end(), points.begin(), cmpByY);
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
    freopen("input_pub.txt", "r", stdin);
    freopen("output.txt", "w", stdout);

    int m;
    cin >> m;

    for (int j = 0; j < m; ++j) {
        cerr << j << endl;

        Polygon p;
        int n, a, b, c;

        cin >> n;
        for (int i = 0; i < n; ++i) {
            cin >> a >> b;
            p.emplace_back(a, b, i);
        }

        auto algo = Kirkpatrick(p);

        cin >> c;
        for (int i = 0; i < c; ++i) {
            cin >> a >> b;
//            cout << (algo.contains({a, b}) ? "INSIDE" : "OUTSIDE") << endl;
        }

    }

    return 0;
}