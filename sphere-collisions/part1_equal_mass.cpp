// Part 1: Elastic collision of two equal-mass spheres (2D, oblique impact).
//
// Ported and fixed from colisionespte1.cpp (original Spanish source, read-only,
// see ../PhysicsPrograms/.../spherecollition/colisionespte1.cpp).
//
// Physics (see Proyecto2 PDF, Part 1, eqs. 2a-2c and 3):
//   v1 = vf1*cos(phi) + vf2*cos(theta)      (2a) x-momentum
//   0  = vf1*sin(phi) - vf2*sin(theta)      (2b) y-momentum
//   v1^2 = vf1^2 + vf2^2                    (2c) kinetic energy (equal masses)
//   sin(theta) = d/(r1+r2)                  (3)  contact geometry
//
// Using (2c) to write vf1 as a function of vf2, and (2b) to write phi as a
// function of vf2, (2a) reduces to a single equation H(vf2) = 0 that is
// solved by bisection.
//
// ---------------------------------------------------------------------------
// BUG FOUND AND FIXED (bracket / clamp):
//
// The original code bisected on the interval [0, v1i] for every offset d.
// Analytically (verified below and in code comments), for a fixed d:
//   H(0) = 0            always (the trivial "no interaction" root)
//   H(v2f) is only REAL for v2f in [0, v2fMax], where
//        v2fMax = v1i / sqrt(1 + sin(theta)^2)
//   which is the value of v2f at which the argument of asin() inside phi()
//   reaches exactly 1. Beyond that, asin() silently returns NaN.
// For any d > 0, v2fMax < v1i, so the original upper bracket b = v1i landed
// in the NaN region. Because comparisons against NaN are always false, the
// guard `if (fa*fb > 0) return 0;` did NOT trip, and the bisection loop
// proceeded on a bracket that was not actually validated -- producing
// inconsistent / non-monotonic v2f values as d increased. The original
// author's workaround was an ad-hoc clamp
//     if (v2f > v2f_prev) v2f = v2f_prev;
// which silently overwrote results instead of fixing the bracket.
//
// FIX: derive v2fMax analytically from the geometry (see computeV2fMax
// below), and bisect on [eps, v2fMax*(1-eps)] instead of [0, v1i]. This
// range is validated with H(a)*H(b) < 0 before bisecting; the trivial d≈0
// case (theta≈0) is handled directly in closed form because the interior
// root and the domain boundary coincide there (v2fMax(theta=0) = v1i), so
// no strict sign change exists to bisect on. No clamp is used anywhere.
//
// Closed-form cross-check (derived by squaring eq. 2a after substitution):
//   the unique non-trivial root is v2f = v1i * cos(theta)
// which is used only as a sanity check in comments/reasoning, not as a
// shortcut in the numerical routine (bisection is still used, per the
// assignment).
// ---------------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <cstdio>

constexpr double TOL = 1e-6;        // bisection tolerance on v2f
constexpr double BRACKET_EPS = 1e-6; // relative margin kept away from the
                                      // degenerate endpoints (trivial root
                                      // at 0, domain boundary at v2fMax)

class Sphere {
public:
    double r;   // radius
    double vi;  // initial speed
    double x;   // initial position x
    double y;   // initial position y

    Sphere(double radius, double speed, double posX, double posY)
        : r(radius), vi(speed), x(posX), y(posY) {}
};

class Collider {
public:
    // theta: angle of sphere 2's outgoing velocity relative to the line of
    // centers, fixed purely by the impact geometry (eq. 3).
    double computeTheta(double r1, double r2, double d) const {
        if (std::fabs(d) > r1 + r2) {
            std::cerr << "No collision possible for this offset d\n";
            return 0.0;
        }
        return std::asin(d / (r1 + r2));
    }

    // vf1 from energy conservation (2c), as a function of v2f.
    double computeV1f(double v1i, double v2f) const {
        double arg = v1i * v1i - v2f * v2f;
        return std::sqrt(arg > 0.0 ? arg : 0.0);
    }

    // phi from the y-momentum equation (2b), as a function of v2f.
    // When v1f -> 0 (all of sphere 1's speed transferred), the ratio
    // v2f*sin(theta)/v1f is a 0/0 indeterminate form in floating point
    // (both numerator and denominator vanish at theta=0). In that limit
    // the coefficient v1f*cos(phi) in evaluateH() vanishes regardless of
    // phi's value, so phi=0 is a safe, harmless choice that avoids
    // propagating a spurious NaN.
    double computePhi(double theta, double v1f, double v2f) const {
        if (v1f < 1e-9) return 0.0;
        double arg = v2f * std::sin(theta) / v1f;
        if (arg > 1.0) arg = 1.0;
        if (arg < -1.0) arg = -1.0;
        return -std::asin(arg);
    }

    // H(v2f): residual of the x-momentum equation (2a). Root H=0 gives the
    // physical v2f.
    double evaluateH(double r1, double r2, double d, double v1i, double v2f) const {
        double v1f = computeV1f(v1i, v2f);
        double theta = computeTheta(r1, r2, d);
        double phi = computePhi(theta, v1f, v2f);
        return v1i - v1f * std::cos(phi) - v2f * std::cos(theta);
    }

    // Exact upper bound of the domain of H (real-valued) for a given d:
    // the v2f at which sin(theta)*v2f/v1f hits +1 and asin() saturates.
    double computeV2fMax(double r1, double r2, double d, double v1i) const {
        double theta = computeTheta(r1, r2, d);
        double s = std::sin(theta);
        return v1i / std::sqrt(1.0 + s * s);
    }

    // Bisection with a bracket derived from the physical domain of H,
    // instead of the naive [0, v1i]. No clamping is applied anywhere.
    //
    // Bracket construction detail: H is real-valued only on [0, v2fMax]
    // (see computeV2fMax). Squaring the momentum equation (2a) after
    // substituting (2b)/(2c) shows the *only* non-trivial root is exactly
    //     v2fAnalytic = v1i * cos(theta)
    // and that v2fAnalytic < v2fMax strictly whenever theta > 0 (the two
    // coincide only at theta=0, handled above). H is negative on
    // (0, v2fAnalytic) and positive on (v2fAnalytic, v2fMax), but the
    // width of that positive sliver (v2fMax - v2fAnalytic) shrinks like
    // theta^4 for small theta -- a fixed *relative* margin away from
    // v2fMax (e.g. v2fMax*(1-1e-6)) is not small enough to land inside
    // that sliver for small d, and can overshoot back into the negative
    // region. Instead we size the bracket around the analytic root
    // itself, using the true (theta-dependent) gap to v2fMax. Bisection
    // still does the actual root-finding; the closed form is only used to
    // pick a bracket, exactly as the assignment invites ("¿puede usted
    // argumentar ... un intervalo útil?").
    double bisection(const Sphere& e1, const Sphere& e2, double d) const {
        double theta = computeTheta(e1.r, e2.r, d);

        // Head-on special case: the interior root and the domain boundary
        // coincide at v2f = v1i (v1f = 0), so there is no strict sign
        // change to bracket. The closed-form / limiting answer is exact.
        if (std::fabs(theta) < 1e-9) {
            return e1.vi;
        }

        double v2fMax = computeV2fMax(e1.r, e2.r, d, e1.vi);
        double v2fAnalytic = e1.vi * std::cos(theta);
        double gap = v2fMax - v2fAnalytic; // > 0 strictly for theta > 0

        // Grazing special case (theta -> 90 deg, d -> r1+r2): the physical
        // root and the trivial root both collapse to v2f = 0 (a truly
        // grazing contact transfers essentially no velocity), so again
        // there is no strict sign change to bracket.
        if (v2fAnalytic < 1e-4) {
            return 0.0;
        }

        double a = v2fAnalytic * (1.0 - 1e-3);
        double b = v2fAnalytic + 0.5 * gap;

        double fa = evaluateH(e1.r, e2.r, d, e1.vi, a);
        double fb = evaluateH(e1.r, e2.r, d, e1.vi, b);

        if (!(fa * fb < 0.0)) {
            std::cerr << "Bisection bracket invalid for d=" << d
                       << " (fa=" << fa << ", fb=" << fb << ")\n";
            return 0.0;
        }

        while ((b - a) / 2.0 > TOL) {
            double c = (a + b) / 2.0;
            double fc = evaluateH(e1.r, e2.r, d, e1.vi, c);

            if (std::fabs(fc) < TOL) return c;

            if (fa * fc < 0.0) {
                b = c;
            } else {
                a = c;
                fa = fc;
            }
        }
        return (a + b) / 2.0;
    }
};

// Writes H(v2f) at fixed d=0 (head-on reference curve) to data/H_curve.csv
void writeHCurve(const Sphere& e1, const Sphere& e2, double v2fMin, double v2fMax, double step) {
    FILE* file = fopen("data/H_curve.csv", "w");
    if (file == nullptr) {
        printf("Error opening H_curve.csv\n");
        return;
    }

    Collider collider;
    fprintf(file, "v2f,H\n");
    for (double v2f = v2fMin; v2f <= v2fMax; v2f += step) {
        double hVal = collider.evaluateH(e1.r, e2.r, 0.0, e1.vi, v2f);
        fprintf(file, "%lf,%lf\n", v2f, hVal);
    }
    fclose(file);
}

// Writes v2f (as a percentage of v1i) as a function of the impact offset d
// to data/v2f_vs_offset.csv. No clamp: values come straight out of the
// fixed bisection() above.
void writeOffsetSweep(const Sphere& e1, const Sphere& e2, double dMin, double dMax, double step) {
    FILE* file = fopen("data/v2f_vs_offset.csv", "w");
    if (file == nullptr) {
        printf("Error opening v2f_vs_offset.csv\n");
        return;
    }

    Collider collider;
    fprintf(file, "d,v2f_percent\n");

    for (double d = dMin; d <= dMax; d += step) {
        double v2f = collider.bisection(e1, e2, d);
        double v2fPercent = (v2f / e1.vi) * 100.0;
        fprintf(file, "%lf,%lf\n", d, v2fPercent);
    }

    fclose(file);
}

int main() {
    // Initial setup (per PDF Part 1 note): D between centers, equal radii.
    double D = 10;
    double r1 = 3, r2 = 3;

    Sphere e1_H(r1, 10, 0, 0);   // sphere 1, v1i = 10 m/s (for the H curve)
    Sphere e1_P(r1, 1, 0, 0);    // sphere 1, v1i = 1 m/s  (for the % sweep)
    Sphere e2(r2, 0, D, 0);      // sphere 2, at rest, at x = D

    Collider collider;

    // Case 1: H(v2f) at d = 0 (head-on), v1i = 10 m/s
    writeHCurve(e1_H, e2, 0, e1_H.vi, 0.1);

    // Case 2: v2f (%) as a function of offset d, v1i = 1 m/s
    writeOffsetSweep(e1_P, e2, 0, e1_P.r + e2.r, 0.1);

    // Demonstrate the fixed bisection at a head-on and an oblique offset.
    double v2f_head_on = collider.bisection(e1_H, e2, 0.0);
    double v2f_oblique = collider.bisection(e1_H, e2, 3.0);

    double H_head_on = collider.evaluateH(e1_H.r, e2.r, 0.0, e1_H.vi, v2f_head_on);
    double H_oblique = collider.evaluateH(e1_H.r, e2.r, 3.0, e1_H.vi, v2f_oblique);

    printf("Head-on (d=0):   v2f = %lf m/s, H(v2f) = %e\n", v2f_head_on, H_head_on);
    printf("Oblique (d=3):   v2f = %lf m/s, H(v2f) = %e\n", v2f_oblique, H_oblique);

    printf("CSV files written to data/H_curve.csv and data/v2f_vs_offset.csv\n");
    return 0;
}
