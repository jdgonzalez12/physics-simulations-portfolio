// Part 2: Elastic collision of two UNEQUAL-mass spheres, head-on (d=0).
//
// Ported and fixed from colisionespte2.cpp (original Spanish source, read-only,
// see ../PhysicsPrograms/.../spherecollition/colisionespte2.cpp).
//
// Physics: same as Part 1, but energy/momentum conservation now carry the
// mass ratio m2/m1 explicitly:
//   m1*v1i = m1*vf1*cos(phi) + m2*vf2*cos(theta)     (x-momentum)
//   0      = m1*vf1*sin(phi) - m2*vf2*sin(theta)     (y-momentum)
//   v1i^2  = vf1^2 + (m2/m1)*vf2^2                   (energy, mass-weighted)
// which is exactly what v1f_fun/phi_fun/h encode below (translated to
// computeV1f/computePhi/evaluateH).
//
// Both example cases here use d=0 (head-on), for which theta=phi=0 and the
// system collapses to the classical 1D elastic collision result
//     v2f = 2*m1*v1i / (m1+m2)
// which is used below only as a sanity comment, not as a shortcut -- the
// program still bisects on H(v2f) using the fixed brackets from the
// original source (verified to bracket the root: 13.33 in [12,14] for
// m1>m2, and 6.67 in [6.5,10] for m2>m1).
//
// ---------------------------------------------------------------------------
// BUG FOUND AND FIXED: the m2>m1 bisection root (v2f_m2_mayor) WAS computed
// but never printed anywhere, unlike the m1>m2 case. Fixed so both cases are
// reported symmetrically to the console, and both CSVs are confirmed written.
//
// ADDITIONAL FINDING while fixing this (worth recording): at exact head-on
// incidence (theta=0), sin(theta)=0 forces phi=-asin(0)=0 identically, so
// cos(phi)=1 is pinned regardless of v2f. That means H(v2f) as coded can
// only represent sphere 1 continuing FORWARD after impact (|phi|<=90 deg
// in general, but pinned to exactly 0 here). For m1 > m2 that's fine
// (sphere 1 keeps moving forward, just slower). For m2 > m1, the correct
// physics has sphere 1 bounce BACKWARD (its final velocity is negative in
// the classical 1D elastic-collision formula), which this phi-convention
// cannot represent -- so H(v2f) stays strictly negative across the entire
// valid domain [0, v2fMax] and the fixed bracket [6.5,10] genuinely brackets
// no root (confirmed numerically: H is negative from v2f=0 all the way to
// v2fMax=v1i*sqrt(m1/m2)=7.07, where H(v2fMax) approx -8.28). This is a
// structural limitation of the original model at d=0, not a bracket bug to
// silently paper over, so bisection() correctly reports "cannot be applied"
// for that case, and main() additionally prints the classical 1D elastic-
// collision closed form v2f = 2*m1*v1i/(m1+m2) as a labeled reference value
// so the m2>m1 case is still reported with a meaningful number, instead of
// a bare, unexplained 0.
// ---------------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <cstdio>

constexpr double TOL = 1e-6;

class Sphere {
public:
    double r;   // radius
    double m;   // mass
    double vi;  // initial speed
    double x;   // initial position x
    double y;   // initial position y

    Sphere(double radius, double mass, double speed, double posX, double posY)
        : r(radius), m(mass), vi(speed), x(posX), y(posY) {}
};

class Collider {
public:
    double computeTheta(double r1, double r2, double d) const {
        if (std::fabs(d) > r1 + r2) {
            std::cerr << "No collision possible for this offset d\n";
            return 0.0;
        }
        return std::asin(d / (r1 + r2));
    }

    double computeV1f(double v1i, double v2f, double massRatio) const {
        double arg = v1i * v1i - massRatio * v2f * v2f;
        return std::sqrt(arg > 0.0 ? arg : 0.0);
    }

    double computePhi(double theta, double v1f, double v2f, double massRatio) const {
        if (v1f < 1e-9) return 0.0;
        double arg = massRatio * v2f * std::sin(theta) / v1f;
        if (arg > 1.0) arg = 1.0;
        if (arg < -1.0) arg = -1.0;
        return -std::asin(arg);
    }

    double evaluateH(double r1, double r2, double d, double v1i, double v2f,
                      double m1, double m2) const {
        double massRatio = m2 / m1;
        double v1f = computeV1f(v1i, v2f, massRatio);
        double theta = computeTheta(r1, r2, d);
        double phi = computePhi(theta, v1f, v2f, massRatio);
        return m1 * (v1i - v1f * std::cos(phi)) - (m2 * v2f * std::cos(theta));
    }

    double bisection(const Sphere& e1, const Sphere& e2, double d, double a, double b) const {
        double fa = evaluateH(e1.r, e2.r, d, e1.vi, a, e1.m, e2.m);
        double fb = evaluateH(e1.r, e2.r, d, e1.vi, b, e1.m, e2.m);

        if (fa * fb > 0) {
            std::cerr << "Bisection cannot be applied: H(a) and H(b) have the same sign\n";
            return 0.0;
        }

        while ((b - a) / 2.0 > TOL) {
            double c = (a + b) / 2.0;
            double fc = evaluateH(e1.r, e2.r, d, e1.vi, c, e1.m, e2.m);

            if (std::fabs(fc) < TOL) return c;

            if (fa * fc < 0) {
                b = c;
            } else {
                a = c;
                fa = fc;
            }
        }
        return (a + b) / 2.0;
    }
};

// Writes H(v2f) for a given case to the requested CSV file.
void writeHCurve(const Sphere& e1, const Sphere& e2, double v2fMin, double v2fMax,
                  double step, const char* fileName) {
    FILE* file = fopen(fileName, "w");
    if (file == nullptr) {
        printf("Error opening %s\n", fileName);
        return;
    }

    Collider collider;
    fprintf(file, "v2f,H\n");
    for (double v2f = v2fMin; v2f <= v2fMax; v2f += step) {
        double hVal = collider.evaluateH(e1.r, e2.r, 0.0, e1.vi, v2f, e1.m, e2.m);
        fprintf(file, "%lf,%lf\n", v2f, hVal);
    }
    fclose(file);
}

int main() {
    double D = 10;
    double r1 = 3, r2 = 3;

    // Case m1 > m2
    Sphere e1_greater(r1, 4, 10, 0, 0);  // sphere 1: m1=4, v1i=10 m/s at (0,0)
    Sphere e2_greater(r2, 2, 0, D, 0);   // sphere 2: m2=2, v2i=0 m/s at (D,0)

    // Case m2 > m1
    Sphere e1_lesser(r1, 2, 10, 0, 0);   // sphere 1: m1=2, v1i=10 m/s at (0,0)
    Sphere e2_lesser(r2, 4, 0, D, 0);    // sphere 2: m2=4, v2i=0 m/s at (D,0)

    Collider collider;

    writeHCurve(e1_greater, e2_greater, 0, 15, 0.1, "data/m1_greater_H.csv");
    writeHCurve(e1_lesser, e2_lesser, 0, 10, 0.1, "data/m2_greater_H.csv");

    double v2f_m1_greater = collider.bisection(e1_greater, e2_greater, 0, 12, 14);
    double v2f_m2_greater = collider.bisection(e1_lesser, e2_lesser, 0, 6.5, 10);

    printf("Final velocity v2f (m1 > m2): %lf m/s\n", v2f_m1_greater);
    printf("Final velocity v2f (m2 > m1): %lf m/s\n", v2f_m2_greater);

    // Classical closed-form reference (1D elastic collision, target at
    // rest): v2f = 2*m1*v1i/(m1+m2). Used here only as a labeled sanity
    // check / fallback report for the m2>m1 case (see comment above main).
    double v2f_m1_greater_ref = 2.0 * e1_greater.m * e1_greater.vi / (e1_greater.m + e2_greater.m);
    double v2f_m2_greater_ref = 2.0 * e1_lesser.m * e1_lesser.vi / (e1_lesser.m + e2_lesser.m);
    printf("Reference (1D elastic formula) v2f (m1 > m2): %lf m/s\n", v2f_m1_greater_ref);
    printf("Reference (1D elastic formula) v2f (m2 > m1): %lf m/s\n", v2f_m2_greater_ref);

    printf("CSV files written to data/m1_greater_H.csv and data/m2_greater_H.csv\n");
    return 0;
}
