// Part 4: Integrative pool-table scenario (see Proyecto2 PDF, Part 4 /
// Figure 6). White ball hits purple ball (swept along y), which then hits
// red ball. Finds the minimum white-ball initial speed so the red ball
// reaches its target position with at least 1 cm/s.
//
// Ported and fixed from colisionespte4.cpp (original Spanish source,
// read-only, see ../PhysicsPrograms/.../spherecollition/colisionespte4.cpp).
// Units, constants (ball mass 250 g, radius parameter 0.06, table 240x120
// cm, friction constant B=0.10) are kept exactly as in the original source
// -- this task's scope is the frictionless/friction split, not re-deriving
// the table geometry.
//
// findRoot() (was encontrar_raiz): a heuristic shrinking-step search for
// the v2f that satisfies the same H(v2f)=0 momentum/energy balance as
// Parts 1-3, applied pairwise to each ball-ball impact along the shot.
// computeFrictionDecay() (was calcular_g): integrates linear drag
// Ff = -B*v over the straight-line segment between two balls and returns
// the surviving fraction of speed (vFinal/vInitial) at the end of the
// segment.
//
// ---------------------------------------------------------------------------
// GAP FOUND AND FIXED: the original program only implemented the
// WITH-FRICTION case (data/graficapte4.csv). The PDF explicitly asks for
// both a frictionless run and a frictional run reported separately. Added
// simulate(bool withFriction) so both data/results_frictionless.csv and
// data/results_with_friction.csv are produced from one shared code path.
// In the frictionless run, computeFrictionDecay() is simply bypassed
// (returns 1.0, i.e. no speed loss along any segment) rather than
// re-deriving a separate physical model.
// ---------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>

constexpr double TOL = 1e-6;
constexpr double B = 0.10;          // friction coefficient (linear drag)
constexpr int MAX_ITER = 10000;

class Ball {
public:
    double m, vi, r, x, y, xf, yf;
    Ball(double mass, double radius, double posX, double posY, double targetX = 0, double targetY = 0)
        : m(mass), vi(0), r(radius), x(posX), y(posY), xf(targetX), yf(targetY) {}
};

double computeTheta(double r1, double r2, double d) {
    if (std::fabs(d) > r1 + r2) return 0;
    return std::asin(d / (r1 + r2));
}

double computeV1f(double v1i, double v2f) {
    double arg = (v1i * v1i) - (v2f * v2f);
    return std::sqrt(arg > 0.0 ? arg : 0.0);
}

double computePhi(double theta, double v1f, double v2f) {
    if (v1f < 1e-9) return 0.0;
    double arg = v2f * std::sin(theta) / v1f;
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;
    return -std::asin(arg);
}

// Heuristic shrinking-step root search for H(v2f)=0 (ported as-is from the
// original encontrar_raiz -- kept close to the source since it is not
// among the identified bugs/gaps for this part).
double findRoot(const Ball& b1, const Ball& b2, double d, double vi) {
    double delta = vi / 100;
    double v2f = delta;
    int iterations = 0;

    while (iterations < MAX_ITER) {
        double theta = computeTheta(b1.r, b2.r, d);
        double v1f = computeV1f(vi, v2f);
        double arg = v1f > 1e-9 ? v2f * std::sin(theta) / v1f : 0.0;

        if ((vi * vi - (v2f * v2f)) < 0 || std::fabs(arg) > 1) {
            v2f -= delta;
            delta /= 2;
        }

        double hVal = b1.m * (vi - v1f * std::cos(computePhi(theta, v1f, v2f))) - (b2.m * v2f * std::cos(theta));
        if (std::fabs(hVal) < TOL) return v2f;

        if (hVal > 0) {
            v2f -= delta;
            delta /= 2;
        }
        v2f += delta;
        iterations++;
    }

    return v2f;
}

// Integrates linear drag Ff=-B*v along the straight segment from (b.x,b.y)
// to (b.xf,b.yf) and returns the surviving fraction of speed vFinal/vInitial.
// When withFriction is false, the segment is frictionless: the ball keeps
// 100% of its speed (fraction = 1.0), and no integration is needed.
double computeFrictionDecay(const Ball& b, bool withFriction) {
    if (!withFriction) return 1.0;

    double d = std::sqrt((b.xf - b.x) * (b.xf - b.x) + (b.yf - b.y) * (b.yf - b.y));
    double tmax = d / b.vi;
    double dt = tmax / 10000;
    double x = 0, vf = b.vi;
    int iterations = 0;

    while (x < d && iterations < MAX_ITER) {
        double xNext = x + (vf * dt) - (0.5 * (B / b.m) * vf * dt * dt);
        if (std::fabs(xNext - x) < TOL) return 0;
        x = xNext;
        vf -= (B / b.m) * vf * dt;
        iterations++;
    }
    return vf / b.vi;
}

// Runs the full white -> purple -> red shot search, sweeping the purple
// ball's y position, and writes (purple_y_cm, white_vi_cm_s) to fileName.
// Returns the best (minimum white speed) position found.
void simulate(bool withFriction, const char* fileName) {
    std::ofstream file(fileName);
    if (!file) {
        std::cerr << "Error opening " << fileName << "\n";
        return;
    }
    file << "purple_y_cm,white_vi_cm_s\n";

    Ball white(0.250, 0.06, 0.40, 0.60);
    Ball purple(0.250, 0.06, 1.20, 0);
    Ball red(0.250, 0.06, 2.00, 1.00, 2.40, 1.20);

    double bestY = 0, bestVi = 1e9;

    for (purple.y = 0.15; purple.y <= 1.20; purple.y += 0.01) {
        double redAngle = std::atan((red.yf - red.y) / (red.xf - red.x));
        purple.xf = red.x - ((purple.r + red.r) * std::cos(redAngle));
        purple.yf = red.y - ((purple.r + red.r) * std::sin(redAngle));
        double purpleAngle = std::atan((purple.yf - purple.y) / (purple.xf - purple.x));
        white.xf = purple.x - ((purple.r + white.r) * std::cos(purpleAngle));
        white.yf = purple.y - ((purple.r + white.r) * std::sin(purpleAngle));
        if (white.yf < 0 || white.xf > 1.20 || white.yf > 1.20) continue;

        white.vi = 0.001;
        int iterations = 0;
        while (iterations < MAX_ITER) {
            double decayWhite = computeFrictionDecay(white, withFriction);
            if (decayWhite == 0) {
                white.vi += 0.001;
                iterations++;
                continue;
            }
            double transferWhitePurple = findRoot(white, purple, std::sin(purpleAngle), 1);
            purple.vi = transferWhitePurple * decayWhite * white.vi;

            double decayPurple = computeFrictionDecay(purple, withFriction);
            if (decayPurple == 0) {
                white.vi += 0.001;
                iterations++;
                continue;
            }
            double transferPurpleRed = findRoot(purple, red, std::sin(redAngle), 1);
            red.vi = decayPurple * transferPurpleRed * transferWhitePurple * decayWhite * white.vi;

            double decayRed = computeFrictionDecay(red, withFriction);
            if (decayRed == 0) {
                white.vi += 0.001;
                iterations++;
                continue;
            }
            if (white.vi * transferWhitePurple * transferPurpleRed * decayWhite * decayPurple * decayRed < 0.01) {
                white.vi += 0.001;
                iterations++;
                continue;
            }
            break;
        }

        double whiteViCms = white.vi * 100;
        if (whiteViCms > 100) whiteViCms = 100;

        file << purple.y * 100 << "," << whiteViCms << "\n";

        if (white.vi < bestVi) {
            bestVi = white.vi;
            bestY = purple.y;
        }
    }

    file.close();
    std::cout << (withFriction ? "[With friction] " : "[Frictionless] ")
              << "Best purple-ball y position: " << bestY * 100 << " cm, "
              << "minimum white-ball speed: " << bestVi * 100 << " cm/s\n";
}

int main() {
    simulate(false, "data/results_frictionless.csv");
    simulate(true, "data/results_with_friction.csv");

    std::cout << "CSV files written to data/results_frictionless.csv and "
                 "data/results_with_friction.csv\n";
    return 0;
}
