#include <tuple>
#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

#define forn(i, n) for(int (i) = 0; (i) < (int)(n); (i)++)
#define forsn(i, s, n) for(int (i) = s; (i) <= (int)(n); (i)++)
#define all(a) (a).begin(), (a).end()

using namespace std;

namespace io {
#include "optimization.h"

#ifndef READ_FILE
#define READ_FILE false
#endif

    void set();
    void writeSpace();
    void writeEndl();

    void set() {
#if READ_FILE
        assert(freopen("input.txt", "r", stdin) != nullptr);
        assert(freopen("output.txt", "w", stdout) != nullptr);
#endif
    }

    void writeSpace() {
        writeChar(' ');
    }

    void writeEndl() {
        writeChar('\n');
    }
}

namespace geom {

    /** DECLARATIONS **/

    inline namespace pt {
        struct point;
        ostream &operator<<(ostream &stream, const point &p);
    }

    /** DEFINITIONS **/

    struct pt::point {
        int x, y;

        point operator-(const point &b) const {
            return {x - b.x, y - b.y};
        }

        point &operator-=(const point &b) {
            x -= b.x;
            y -= b.y;
            return *this;
        }

        long long operator*(const point &b) const {
            return this->x * 1ll * b.x + this->y * 1ll * b.y;
        }

        long long operator%(const point &b) const {
            return this->x * 1ll * b.y - this->y * 1ll * b.x;
        }

        bool operator<(const point &b) const {
            return this->x < b.x || (this->x == b.x && this->y < b.y);
        }

        bool operator==(const point &b) const {
            return this->x == b.x && this->y == b.y;
        }

        double dist() const {
            return sqrt(pow((double) x, 2) + pow((double) y, 2));
        }

        double atan() const {
            return atan2((double) y, (double) (x == 0 && y == 0) ? 1 : x);
        }
    };

    ostream &pt::operator<<(ostream &stream, const point &p) {
        stream << "(" << p.x << ", " << p.y << ")";
        return stream;
    }

}

namespace qst {
    template<class T>
    class AddMaxSegmTree {
        int n;
        vector<T> s;
        vector<pair<T, int>> m;

        tuple<int, int, int> nextNodes(int v, int lt, int rt) {
            int mt = (lt + rt) >> 1;
            int lv = v << 1;
            int rv = lv + 1;
            return make_tuple(mt, lv, rv);
        };

        void build(const vector<T> &a, int v, int lt, int rt) {
            if (lt == rt) {
                s[v] = a[lt];
                m[v] = make_pair(a[lt], lt);
            } else {
                int mt, lv, rv;
                tie(mt, lv, rv) = nextNodes(v, lt, rt);
                build(a, lv, lt, mt);
                build(a, rv, mt + 1, rt);
                m[v] = std::max(m[lv], m[rv]);
            }
        }

        void add(int v, int lt, int rt, int l, int r, T delta) {
            if (l > r) {
                return;
            }
            if (lt == l && rt == r) {
                s[v] += delta;
                m[v].first += delta;
                return;
            }
            int mt, lv, rv;
            tie(mt, lv, rv) = nextNodes(v, lt, rt);
            add(lv, lt, mt, l, min(r, mt), delta);
            add(rv, mt + 1, rt, std::max(l, mt + 1), r, delta);
            m[v] = std::max(m[lv], m[rv]);
            m[v].first += s[v];
        }

        pair<T, int> max(int v, int lt, int rt, int l, int r) {
            if (l > r) {
                return make_pair(numeric_limits<T>::min(),
                                 numeric_limits<int>::min());
            }
            if (lt == l && rt == r) {
                return m[v];
            }
            int mt, lv, rv;
            tie(mt, lv, rv) = nextNodes(v, lt, rt);
            auto ans = std::max(max(lv, lt, mt, l, min(r, mt)),
                                max(rv, mt + 1, rt, std::max(l, mt + 1), r));
            ans.first += s[v];
            return ans;
        }

    public:
        explicit AddMaxSegmTree(const vector<T> &a) {
            n = (int) a.size();
            s.resize((size_t) n << 2, 0);
            m.resize((size_t) n << 2);
            build(a, 1, 0, n - 1);
        }

        void add(int l, int r, T delta) {
            add(1, 0, n - 1, l, r, delta);
        }

        pair<T, int> max(int l, int r) {
            return max(1, 0, n - 1, l, r);
        }
    };
}

void solve();

int main() {
    io::set();
    solve();
    io::flush();
}

void solve() {
    using namespace geom::pt;
    using qst::AddMaxSegmTree;
    int s = io::readInt(), w = io::readInt();
    int n = io::readInt();
    vector<point> r;
    r.reserve(static_cast<unsigned long>(n));
    forn(_, n) {
        r.push_back({io::readInt(), io::readInt()});
    }
    int ss[4] = {-s, -s, s, s};
    int ws[4] = {-w, w, -w, w};
    auto shift = (int) (4 * (1e+4));
    int ans_num = 0;
    point ans_p{0, 0};
    forn(k, 4) {
        s = ss[k];
        w = ws[k];
        vector<tuple<int, int, int, int>> a;
        a.reserve(static_cast<unsigned long>(2 * n)); // NOLINT
        forn(i, n) {
            auto p = r[i];
            a.emplace_back(p.y - w, 0, p.x - s, p.x);
            a.emplace_back(p.y, 1, p.x - s, p.x);
        }
        sort(all(a));
        auto xq = AddMaxSegmTree<int>(vector<int>((size_t) (2 * shift + 1), 0));
        forn(i, a.size()) {
            int y, c, lx, rx;
            tie(y, c, lx, rx) = a[i];
            if (c == 0) {
                xq.add(lx + shift, rx + shift, 1);
                auto m = xq.max(0, 2 * shift);
                if (m.first > ans_num) {
                    ans_num = m.first;
                    ans_p = point{m.second - shift + s, y + w};
                }
            } else {
                xq.add(lx + shift, rx + shift, -1);
            }
        }
    }
    io::writeInt(ans_num);
    io::writeEndl();
    io::writeInt(ans_p.x);
    io::writeSpace();
    io::writeInt(ans_p.y);
    io::writeEndl();
}