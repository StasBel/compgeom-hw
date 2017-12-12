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
    class AddOneSumSegmTree {
        int n;
        vector<T> t;

        void build(const vector<T> &a, int v, int lt, int rt) {
            if (lt == rt) {
                t[v] = a[lt];
            } else {
                int mt = (lt + rt) >> 1;
                int lv = v << 1;
                int rv = lv + 1;
                build(a, lv, lt, mt);
                build(a, rv, mt + 1, rt);
                t[v] = t[lv] + t[rv];
            }
        }

        void add(int v, int lt, int rt, int pos, T delta) {
            if (lt == rt) {
                t[v] += delta;
            } else {
                int mt = (lt + rt) >> 1;
                int lv = v << 1;
                int rv = lv + 1;
                if (pos <= mt) {
                    add(lv, lt, mt, pos, delta);
                } else {
                    add(rv, mt + 1, rt, pos, delta);
                }
                t[v] = t[lv] + t[rv];
            }
        }

        T sum(int v, int lt, int rt, int l, int r) {
            if (l > r) {
                return T(0);
            }
            if (lt == l && rt == r) {
                return t[v];
            }
            int mt = (lt + rt) >> 1;
            int lv = v << 1;
            int rv = lv + 1;
            return sum(lv, lt, mt, l, min(r, mt))
                   + sum(rv, mt + 1, rt, max(l, mt + 1), r);
        }

    public:
        explicit AddOneSumSegmTree(const vector<T> &a) {
            n = (int) a.size();
            t.resize((size_t) n << 2);
            build(a, 1, 0, n - 1);
        }

        void add(int pos, T delta) {
            add(1, 0, n - 1, pos, delta);
        }

        T sum(int l, int r) {
            return sum(1, 0, n - 1, l, r);
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
    using qst::AddOneSumSegmTree;

    int n;
    while ((n = io::readInt())) {
        auto xq = AddOneSumSegmTree<int>(vector<int>(32001, 0));
        auto levels = vector<int>((size_t) n, 0);
        forn(_, n) {
            auto p = point{io::readInt(), io::readInt()};
            levels[xq.sum(0, p.x)] += 1;
            xq.add(p.x, 1);
        }
        forn(i, n) {
            io::writeInt(levels[i]);
            io::writeEndl();
        }
    }
}