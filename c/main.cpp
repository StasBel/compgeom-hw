#include "optimization.h"
#include <iostream>
#include <vector>
#include <cmath>

#define forn(i, n) for(int (i) = 0; (i) < (int)(n); (i)++)
#define forsn(i, s, n) for(int (i) = s; (i) <= (int)(n); (i)++)
#define all(a) (a).begin(), (a).end()

using namespace std;

#define READ_FILE true

void set_io() {
    assert(freopen("input.txt", "r", stdin) != nullptr);
    assert(freopen("output.txt", "w", stdout) != nullptr);
}

inline void writeEndl() {
    writeChar('\n');
}

struct point {
    int x, y;

    inline point operator-(const point &b) const {
        return {x - b.x, y - b.y};
    }

    inline point &operator-=(const point &b) {
        x -= b.x;
        y -= b.y;
        return *this;
    }

    inline long operator*(const point &b) const {
        return this->x * 1ll * b.x + this->y * 1ll * b.y;
    }

    inline long long operator%(const point &b) const {
        return this->x * 1ll * b.y - this->y * 1ll * b.x;
    }

    inline bool operator<(const point &b) const {
        return this->x < b.x || (this->x == b.x && this->y < b.y);
    }

    inline double dist() const {
        return sqrt(pow(x, 2) + pow(y, 2));
    }

    inline double atan() const {
        return atan2(y, (x == 0 && y == 0) ? 1 : x);
    }
};

std::ostream &operator<<(std::ostream &stream, const point &p) {
    stream << "(" << p.x << ", " << p.y << ")";
    return stream;
}

/**
 * @return Return 1 for right turn, 0 for on line, -1 for left turn.
 */
inline int turn(point p, point q, point r) {
    auto vec = (q - p) % (r - p);
    return (vec < 0.) - (vec > 0.);
}

inline vector<point> graham_ch(vector<point> P) {
    auto n = (int) P.size();
    if (n <= 1) {
        return P;
    }
    auto it = min_element(all(P), [](const point &p, const point &q) {
        return p.y < q.y || (p.y == q.y && p.x > q.x);
    });
    iter_swap(it, P.begin());
    const point &pivot = P[0];
    sort(++all(P), [pivot](const point &p, const point &q) {
        auto pp = p - pivot, qp = q - pivot;
        auto pa = pp.atan(), qa = qp.atan();
        return pa < qa || (pa == qa && pp.dist() < qp.dist());
    });
    int k = 1;
    forsn(i, 2, n - 1) {
        while (k > 0 && turn(P[k - 1], P[k], P[i]) >= 0) {
            k--;
        }
        swap(P[i], P[k + 1]);
        k++;
    }
    P.resize((size_t) min(k + 1, n));
    return P;
}

void solve() {
    int n = readInt();
    writeInt(n);
}

int main() {
#if READ_FILE
    set_io();
#endif
    solve();
#if READ_FILE
    flush();
#endif
}