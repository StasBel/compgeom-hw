#include "optimization.h"
#include <iostream>
#include <vector>
#include <cmath>

#define forn(i, n) for(int (i) = 0; (i) < (int)(n); (i)++)
#define forsn(i, s, n) for(int (i) = s; (i) <= (int)(n); (i)++)
#define all(a) (a).begin(), (a).end()

using namespace std;

#define READ_FILE false

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
        return sqrt(pow((double) x, 2) + pow((double) y, 2));
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

inline double turn_len(const point &p, const point &q, const point &r, double l) {
    auto a = q - r, b = q - p;
    auto dp = (double) (a * b);
    dp /= a.dist();
    dp /= b.dist();
    auto angle = M_PI - acos(dp);
    return angle * l;
}

void solve() {
    int n = readInt();
    double l = readDouble();
    vector<point> P;
    P.reserve((size_t) n);
    forn(_, n) {
        P.push_back({readInt(), readInt()});
    }
    auto Q = graham_ch(P);
    auto m = Q.size();
    double len = 0;
    forn(i, m - 2) {
        len += turn_len(Q[i], Q[i + 1], Q[i + 2], l);
    }
    forn(i, m - 1) {
        len += (Q[i + 1] - Q[i]).dist();
    }
    len += turn_len(Q[m - 2], Q[m - 1], Q[0], l);
    len += turn_len(Q[m - 1], Q[0], Q[1], l);
    len += (Q[0] - Q[m - 1]).dist();
    writeDouble(len, 9);
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