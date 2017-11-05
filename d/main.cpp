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
    inline namespace op {
        int turn(const point &p, const point &q, const point &r);
        double angle(const point &p, const point &q, const point &r);
        double area(const vector<point> &P);
    }
    namespace in {
        class Pangle;
    }
    namespace ch {
        bool graham(vector<point> &P);
        int graham_quick(vector<point> &P);
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

    /**
     * @return Return 1 for right turn, 0 for on line, -1 for left turn.
     */
    int op::turn(const point &p, const point &q, const point &r) {
        auto vec = (q - p) % (r - p);
        return (vec < 0.) - (vec > 0.);
    }

    double op::area(const vector<point> &P) {
        double answer = 0;
        forn(i, P.size()) {
            auto p1 = (i ? P[i - 1] : P.back()), p2 = P[i];
            answer += (double) ((p1.x - p2.x) * (p1.y + p2.y));
        }
        return abs(answer) / 2.;
    }

    /**
     * @return Angle between p-q and p-r, anti-cloakwise;
     */
    double op::angle(const point &p, const point &q, const point &r) {
        auto a = q - p, b = r - p;
        auto dp = (double) (a * b);
        dp /= a.dist();
        dp /= b.dist();
        return acos(dp);
    }

    class in::Pangle {
    private:
        vector<point> P;
        point pivot{};
        vector<double> ang;
    public:
        explicit Pangle(const vector<point> &_P) {
            P = _P;
            auto pivot_it = min_element(all(P));
            pivot = *pivot_it;
            rotate(P.begin(), pivot_it, P.end());
            P.erase(P.begin());
            forn(i, P.size()) {
                P[i] -= pivot;
            }
            ang.reserve(P.size());
            forn(i, P.size()) {
                ang.push_back(P[i].atan());
            }
            if (ang[0] > ang[P.size() - 1]) {
                reverse(all(P));
                reverse(all(ang));
            }
        }

        int check(const point &_q) {
            auto q = _q;
            q -= pivot;

            int answer = -1;
            if (q.x >= 0) {
                if (q.x == 0 && q.y == 0) {
                    answer = 0;
                } else {
                    auto angle = q.atan();
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

            return answer;
        }
    };

    bool ch::graham(vector<point> &P) {
        auto n = (int) P.size();

        if (n <= 1) {
            return false;
        }

        auto it = min_element(all(P), [](const point &p, const point &q) {
            return p.y < q.y || (p.y == q.y && p.x > q.x);
        });
        iter_swap(it, P.begin());

        auto pivot = P[0];
        vector<pair<pair<double, double>, point> > Q;
        transform(all(P), back_inserter(Q), [&pivot](const point &p) {
            auto pp = p - pivot;
            return make_pair(make_pair(pp.atan(), pp.dist()), p);
        });
        sort(++all(Q));
        P.clear();
        transform(all(Q), back_inserter(P), [](const pair<pair<double, double>, point> &p) {
            return p.second;
        });

        int k = 1;
        vector<point> I, B;
        forsn(i, 2, n - 1) {
            while (k > 0 && turn(P[k - 1], P[k], P[i]) > 0) {
                while (!B.empty() && turn(P[k - 1], P[k], B.back()) == 0) {
                    I.emplace_back(B.back());
                    B.pop_back();
                }
                I.push_back(P[k]);
                k--;
            }
            while (k > 0 && turn(P[k - 1], P[k], P[i]) == 0) {
                B.push_back(P[k]);
                k--;
            }
            swap(P[i], P[k + 1]);
            k++;
        }
        while (!I.empty() && turn(P[0], P[k], I.back()) == 0) {
            I.pop_back();
        }
        int m = min(k + 1, n);
        bool not_v = turn(P[0], P[1], P[m - 1]) != 0;
        P = I;
        return not_v;
    }

    int ch::graham_quick(vector<point> &P) {
        auto n = (int) P.size();

        if (n <= 1) {
            return 1;
        }

        auto it = min_element(all(P), [](const point &p, const point &q) {
            return p.y < q.y || (p.y == q.y && p.x > q.x);
        });
        iter_swap(it, P.begin());

        auto pivot = P[0];
        vector<pair<pair<double, double>, point> > Q;
        transform(all(P), back_inserter(Q), [&pivot](const point &p) {
            auto pp = p - pivot;
            return make_pair(make_pair(pp.atan(), pp.dist()), p);
        });
        sort(++all(Q));
        P.clear();
        transform(all(Q), back_inserter(P), [](const pair<pair<double, double>, point> &p) {
            return p.second;
        });

        int k = 1;
        forsn(i, 2, n - 1) {
            while (k > 0 && turn(P[k - 1], P[k], P[i]) >= 0) {
                k--;
            }
            swap(P[i], P[k + 1]);
            k++;
        }
        return min(k + 1, n);
    }
}

void solve();

int main() {
    io::set();
    solve();
    io::flush();
}

void solve() {
    auto n = io::readInt();
    using namespace geom::pt;
    vector<point> P;
    P.reserve((size_t) n);
    forn(_, n) {
        P.push_back({io::readInt(), io::readInt()});
    }

    int c = 0;
    using geom::ch::graham;
    using geom::op::turn;
    while (P.size() >= 3) {
        if (graham(P)) {
            c += 1;
        }
    }
    io::writeInt(c);
}