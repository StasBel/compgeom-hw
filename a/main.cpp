#include "optimization.h"
#include <iostream>
#include <vector>
#include <cmath>

#define forn(i, n) for(int (i) = 0; (i) < (int)(n); (i)++)
#define forsn(i, s, n) for(int (i) = s; (i) <= (int)(n); (i)++)
#define all(a) (a).begin(), (a).end()

using namespace std;

void set_io() {
    assert(freopen("input.txt", "r", stdin) != nullptr);
    assert(freopen("output.txt", "w", stdout) != nullptr);
}

struct point {
    int x, y;

    point operator-(const point &b) const {
        return {this->x - b.x, this->y - b.y};
    }

    point &operator-=(const point &b) {
        this->x -= b.x;
        this->y -= b.y;
        return *this;
    }

    long operator*(const point &b) const {
        return this->x * 1ll * b.x + this->y * 1ll * b.y;
    }

    long long operator%(const point &b) const {
        return this->x * 1ll * b.y - this->y * 1ll * b.x;
    }

    bool operator<(const point &b) const {
        return this->x < b.x || (this->x == b.x && this->y < b.y);
    }
};

std::ostream &operator<<(std::ostream &stream, const point &p) {
    stream << "{" << p.x << ", " << p.y << "}";
    return stream;
}

/**
 * @return 1 for right, 0 for on-line, -1 for left
 */
inline int turn(point p, point q, point r) {
    auto vec = (q - p) % (r - p);
    return (vec < 0.) - (vec > 0.);
}

long double atan3(const point &p) {
    return atan2((long double) p.y, (long double) ((p.y == 0 && p.x == 0) ? 1 : p.x));
}

void solve() {
    // read
    int n = readInt();
    vector<point> P;
    P.reserve((size_t) n);
    forn(_, n) {
        P.push_back({readInt(), readInt()});
    }

    // preprocess
    auto it_pivot = min_element(all(P));
    point pivot = *it_pivot;
    rotate(P.begin(), it_pivot, P.end());
    P.erase(P.begin());
    forn(i, P.size()) {
        P[i] -= pivot;
    }
    vector<long double> ang;
    ang.reserve(P.size());
    forn(i, P.size()) {
        ang.push_back(atan3(P[i]));
    }
    if (ang[0] > ang[P.size() - 1]) {
        reverse(all(P));
        reverse(all(ang));
    }

    // queries
    int k = readInt();
    forn(_, k) {
        point q = {readInt(), readInt()};
        q -= pivot;

        int answer = -1;
        if (q.x >= 0) {
            if (q.x == 0 && q.y == 0) {
                answer = 0;
            } else {
                long double angle = atan3(q);
                if ((angle >= ang[0]) && (angle <= ang[ang.size() - 1])) {
                    auto it_up = upper_bound(all(ang), angle);
                    point r{0, 0};
                    if (it_up != ang.end()) {
                        r = P[it_up - ang.begin()];
                    }
                    point p = P[it_up - ang.begin() - 1];
                    int turn1 = turn(p, q, r);

                    auto it_lo = lower_bound(all(ang), angle);
                    if (it_lo == ang.begin()) {
                        p = {0, 0};
                    } else {
                        p = P[it_lo - ang.begin() - 1];
                    }
                    r = P[it_lo - ang.begin()];
                    int turn2 = turn(p, q, r);
                    int sum_turn = turn1 + turn2;
                    if (sum_turn == 2) {
                        answer = 1;
                    } else if (turn1 == -1 || turn2 == -1) {
                        answer = -1;
                    } else {
                        answer = 0;
                    }
                }
            }
        }

        switch (answer) {
            case -1:
                writeWord("OUTSIDE\n");
                break;
            case 0:
                writeWord("BORDER\n");
                break;
            case 1:
                writeWord("INSIDE\n");
                break;
            default:
                return;
        }
    }
    flush();
}

int main() {
    // set_io();
    solve();
}