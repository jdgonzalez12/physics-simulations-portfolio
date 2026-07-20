// Part 3: Elastic collision of two unequal masses, oblique incidence angle
// alpha (see Proyecto2 PDF, Part 3 and Figure 4/5).
//
// Ported and fixed from colisionespte3.cpp (original Spanish source, read-only,
// see ../PhysicsPrograms/.../spherecollition/colisionespte3.cpp).
//
// Model (unchanged from the original source): sphere 1's initial velocity is
// still along the x-axis, but sphere 2's initial position is rotated by
// alpha away from the x-axis, at fixed distance D from sphere 1:
//     Q = (D*cos(alpha), D*sin(alpha))
// The x-momentum equation (2a) is modified so the along-line-of-centers
// component of sphere 1's incoming momentum is reduced by cos(alpha):
//     m1*(v1i*cos(alpha) - vf1*cos(phi)) = m2*vf2*cos(theta)
// exactly as in evaluateH() below. As in the original source, d=0 is used
// in the theta()/H() calls for every alpha (this was true even for the
// single alpha=30deg case in the original code, not something introduced
// here) -- i.e. this models a straight-line, center-aligned approach whose
// effective closing speed is attenuated by cos(alpha), consistent with a
// head-on impact happening along the (rotated) P-Q line. With d=0, theta
// and phi both vanish identically (see Part 1/2). Unlike Part 2 though,
// v0=v1i*cos(alpha) only enters the momentum term here -- energy
// conservation still uses the true v1i (alpha attenuates how much of
// sphere 1's momentum acts along the line of centers, not its actual
// kinetic energy) -- so the reduced equation is a genuine quadratic in
// v2f, not Part 2's simple linear relation. Its positive-branch root
// (derived in bisection() below) is used only to build a numerically safe
// bisection bracket per alpha value, the same technique used to fix
// Part 1's bracket bug -- bisection still does the actual root-finding.
//
// ---------------------------------------------------------------------------
// GAP FOUND AND FIXED: the original program only ever evaluated one
// hardcoded alpha = 30 deg (pi/6). Per the PDF (Part 3 note + Figure 5),
// the valid range of alpha is bounded by the grazing-impact geometry: as
// alpha grows, sphere 2's position Q=(D cos a, D sin a) moves away from
// sphere 1's straight-line path (the x-axis), and the perpendicular
// distance from Q to that path is D*sin(alpha). A collision is only
// geometrically possible while that perpendicular distance does not
// exceed r1+r2 (exactly the same grazing condition as eq. 3, sin(theta) =
// d/(r1+r2), now applied to the *initial* geometry rather than a variable
// offset). This gives:
//     alpha_max = asin((r1+r2) / D)
// which is precisely the isosceles-triangle relation illustrated in
// Figure 5 (legs of length D from P to the two extreme tangent
// directions, half-angle alpha_max, opposite side r1+r2 on each half).
// With D=10 and r1=r2=3, alpha_max = asin(0.6) = 36.87 deg -- notably
// *larger* than the original hardcoded 30 deg, confirming that value was
// an arbitrary interior sample rather than a derived bound.
//
// The program below sweeps alpha over [0, alpha_max) (grazing alpha_max
// itself is excluded, same degenerate-bracket reasoning as Part 1's d ->
// r1+r2 edge) with 80 steps, for both mass-ratio cases, writing
// data/m1_greater_alpha.csv and data/m2_greater_alpha.csv with columns
// alpha_deg, alpha_rad, v2f.
//
// Counter-intuitive but verified result for m1>m2: v2f INCREASES slightly
// with alpha (13.33 m/s at alpha=0 up to 14.07 m/s at alpha_max, versus a
// fixed energy ceiling v2fMax=v1i/sqrt(massRatio)=14.14 m/s). This is
// because alpha only weakens the MOMENTUM constraint while sphere 1's
// total kinetic energy budget stays fixed at v1i -- a weaker momentum
// constraint against a fixed energy budget lets the root sit closer to
// the energy ceiling, not further from it.
//
// NOTE (carried over from the Part 2 finding): since theta=phi=0 for every
// alpha here (d=0 fixed, see above), the same structural limitation found
// in Part 2 applies for *every* alpha whenever m2 > m1, not just alpha=0.
// Solving the quadratic and checking its sign against the ORIGINAL
// (unsquared) relation v1f = v0 - massRatio*v2f shows v1f would have to be
// negative (sphere 1 bouncing backward) at both alpha=0 and alpha near
// alpha_max, and by continuity everywhere in between -- confirmed
// numerically across the full sweep. This phi-convention cannot represent
// a backward bounce, so H(v2f) has no interior root anywhere in the
// m2>m1 sweep. bisection() detects this up front (v1fSigned < 0 check)
// and returns 0.0 immediately, without spamming 80 identical
// bracket-invalid warnings -- the m2>m1 CSV is therefore a flat v2f=0
// curve, an honest reflection of the model's limitation rather than a bug.
// ---------------------------------------------------------------------------

#include <iostream>
#include <cmath>
#include <cstdio>

constexpr double TOL = 1e-6;
constexpr double PI = 3.14159265358979323846;

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

    // H(v2f) with the oblique-incidence correction v1i*cos(alpha).
    double evaluateH(double r1, double r2, double d, double v1i, double v2f,
                      double m1, double m2, double alpha) const {
        double massRatio = m2 / m1;
        double v1f = computeV1f(v1i, v2f, massRatio);
        double theta = computeTheta(r1, r2, d);
        double phi = computePhi(theta, v1f, v2f, massRatio);
        return m1 * (v1i * std::cos(alpha) - v1f * std::cos(phi)) - (m2 * v2f * std::cos(theta));
    }

    // Bisection with a bracket sized from a closed-form estimate of the
    // root of the theta=0 (d=0) reduction, the same strategy used to fix
    // Part 1's bracket bug: a fixed bracket like [12,14] only works for
    // the single alpha it was tuned to and breaks as alpha sweeps.
    //
    // Careful derivation: at theta=phi=0, H(v2f)=0 means
    //     m1*(v0 - v1f) = m2*v2f,   v1f = sqrt(v1i^2 - massRatio*v2f^2)
    // Note v0 = v1i*cos(alpha) appears ONLY in the momentum term; v1f's
    // energy relation still uses the true v1i (alpha does not change
    // sphere 1's actual kinetic energy, only the fraction of its momentum
    // directed along the line of centers). Squaring the signed relation
    // v1f = v0 - massRatio*v2f and substituting gives a genuine quadratic
    // in v2f (NOT the simple linear result of Part 2, which only held
    // because v0=v1i there):
    //     massRatio*(massRatio+1)*v2f^2 - 2*v0*massRatio*v2f + (v0^2-v1i^2) = 0
    // whose physically relevant root (positive v2f, positive v1f) is the
    // "+" branch used below.
    double bisection(const Sphere& e1, const Sphere& e2, double d, double alpha) const {
        double massRatio = e2.m / e1.m;
        double v0 = e1.vi * std::cos(alpha);
        double v1i = e1.vi;

        double discriminant = ((massRatio + 1.0) * v1i * v1i - v0 * v0) / massRatio;
        double v2fAnalytic = (v0 + std::sqrt(discriminant > 0 ? discriminant : 0.0)) / (massRatio + 1.0);

        // Domain boundary where v1f -> 0: v1f^2 = v1i^2 - massRatio*v2f^2
        // depends on the TRUE v1i, not on the alpha-attenuated v0.
        double v2fMax = v1i / std::sqrt(massRatio);

        // Reject roots that don't satisfy the ORIGINAL signed relation
        // v1f = v0 - massRatio*v2f (squaring can introduce a spurious
        // root where this would require v1f < 0, i.e. sphere 1 bouncing
        // backward -- unrepresentable at theta=0, same finding as Part 2).
        double v1fSigned = v0 - massRatio * v2fAnalytic;

        if (v2fAnalytic < 1e-4 || v2fMax - v2fAnalytic < 1e-9 || v1fSigned < 0.0) {
            return 0.0;
        }

        double gap = v2fMax - v2fAnalytic;
        double a = v2fAnalytic * (1.0 - 1e-3);
        double b = v2fAnalytic + 0.5 * gap;

        double fa = evaluateH(e1.r, e2.r, d, e1.vi, a, e1.m, e2.m, alpha);
        double fb = evaluateH(e1.r, e2.r, d, e1.vi, b, e1.m, e2.m, alpha);

        if (!(fa * fb < 0.0)) {
            std::cerr << "Bisection bracket invalid for alpha=" << alpha
                       << " (fa=" << fa << ", fb=" << fb << ")\n";
            return 0.0;
        }

        while ((b - a) / 2.0 > TOL) {
            double c = (a + b) / 2.0;
            double fc = evaluateH(e1.r, e2.r, d, e1.vi, c, e1.m, e2.m, alpha);

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

// Sweeps alpha in [0, alphaMax) and writes (alpha_deg, alpha_rad, v2f) to
// the requested CSV file. For m1 >= m2, v2f comes from genuine bisection on
// H(v2f), which finds a valid root at every alpha in the sweep (verified).
// For m2 > m1, the theta=0 model has no interior root at ANY alpha (proven
// in the header note: the signed v1f = v0 - massRatio*v2f is negative at
// both alpha=0 and alpha->alphaMax, and by continuity in between too), so
// bisection() correctly and consistently returns 0.0 there without any
// bracket-invalid spam (the rejection happens up front via the v1fSigned
// check, before any bisection iterations are attempted).
void writeAlphaSweep(const Sphere& e1, const Sphere& e2, double alphaMax,
                      int steps, const char* fileName) {
    FILE* file = fopen(fileName, "w");
    if (file == nullptr) {
        printf("Error opening %s\n", fileName);
        return;
    }

    Collider collider;
    fprintf(file, "alpha_deg,alpha_rad,v2f\n");

    for (int i = 0; i < steps; ++i) {
        double alpha = alphaMax * (static_cast<double>(i) / steps); // [0, alphaMax)
        double v2f = collider.bisection(e1, e2, 0.0, alpha);
        double alphaDeg = alpha * 180.0 / PI;
        fprintf(file, "%lf,%lf,%lf\n", alphaDeg, alpha, v2f);
    }

    fclose(file);
}

int main() {
    double D = 10;
    double r1 = 3, r2 = 3;

    // Grazing-impact bound derived from the isosceles-triangle geometry of
    // Figure 5: sin(alpha_max) = (r1+r2)/D.
    double alphaMax = std::asin((r1 + r2) / D);
    printf("Valid alpha range: [0, %lf] rad = [0, %lf] deg\n",
           alphaMax, alphaMax * 180.0 / PI);

    // Case m1 > m2. Sphere 2's position is only used for documentation
    // here (evaluateH uses d=0 for every alpha, see header comment); it is
    // still placed at the rotated reference angle for consistency with the
    // original source's geometric setup.
    Sphere e1_greater(r1, 4, 10, 0, 0);
    Sphere e2_greater(r2, 2, 0, D * std::cos(alphaMax / 2), D * std::sin(alphaMax / 2));

    // Case m2 > m1
    Sphere e1_lesser(r1, 2, 10, 0, 0);
    Sphere e2_lesser(r2, 4, 0, D * std::cos(alphaMax / 2), D * std::sin(alphaMax / 2));

    const int STEPS = 80;
    writeAlphaSweep(e1_greater, e2_greater, alphaMax, STEPS, "data/m1_greater_alpha.csv");
    writeAlphaSweep(e1_lesser, e2_lesser, alphaMax, STEPS, "data/m2_greater_alpha.csv");

    Collider collider;
    double alphaSample = PI / 6.0; // 30 deg, the original single sample point
    double v2f_m1_greater = collider.bisection(e1_greater, e2_greater, 0.0, alphaSample);
    double v2f_m2_greater = collider.bisection(e1_lesser, e2_lesser, 0.0, alphaSample);

    printf("At alpha = 30 deg: v2f (m1 > m2) = %lf m/s (via bisection)\n", v2f_m1_greater);
    printf("At alpha = 30 deg: v2f (m2 > m1) = %lf m/s (no interior root under this model, see header note)\n", v2f_m2_greater);

    printf("CSV files written to data/m1_greater_alpha.csv and data/m2_greater_alpha.csv\n");
    return 0;
}
