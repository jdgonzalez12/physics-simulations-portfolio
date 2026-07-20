// part2_velocity_acceleration.cpp
//
// Part 2: Numerical velocity and acceleration estimation for points A and B
// of the 5-bar oscillating linkage, via central finite differences applied
// to the position time series produced the same way as in Part 1.
//
// theta2 is swept over two full periods [0, 4*pi] with a constant step. This
// implicitly assumes constant angular velocity omega = 1 rad/s, so theta2
// doubles as a dimensionless time proxy (t = theta2 / omega = theta2 when
// omega = 1 rad/s) -- see D_THETA2 below and the notebook axis labels.
//
// --- BUG FIX (relative to the original barras2.cpp) -----------------------
// The original code defined a finite-difference step constant
//   #define DTH2 1e-6
// that was completely disconnected from the actual angular increment used by
// the sweep loop, (4*pi)/N with N=1000 (~= 0.012566). Central differences
// divided by 2*DTH2 (~2e-6) instead of 2*(4*pi/N) (~0.025133) inflated the
// computed velocities by a factor of roughly (4*pi/N)/DTH2 ~= 1.26e4. The
// acceleration formula divides by DTH2*DTH2, so the same mismatch is squared,
// inflating accelerations by roughly [(4*pi/N)/DTH2]^2 ~= 1.6e8. This is
// exactly why the original program's velocity/acceleration output landed
// around 1e6-1e9 instead of the O(10-500) range shown in the assignment's
// reference Fig. 3.
//
// Fix: derive a single constant D_THETA2 = (4*pi)/N from the same N used by
// the sweep loop, and reuse it as both the loop increment and the finite-
// difference step, so there is one source of truth instead of two constants
// (DTH2 and the loop step) that had to (but didn't) stay in sync.
// ---------------------------------------------------------------------------

#define _USE_MATH_DEFINES // needed for M_PI on MinGW/MSVC
#include <iostream>
#include <cstdio>
#include <cmath>

#define TOL 1e-6
#define MAX_ITER 1000

const int N = 1000;
// Angular step used by the theta2 sweep below, over two periods [0, 4*pi].
// Reused as the finite-difference step so velocity/acceleration formulas are
// consistent with the actual sampling interval (single source of truth).
// Assumes constant angular velocity omega = 1 rad/s, so theta2 acts as time.
const double D_THETA2 = (4.0 * M_PI) / N;

class Mechanism {
public:
    double r1, r2, r3, r4, r5, alpha;
    double theta2, theta3, theta4;
    double k1, k2, k3;

    Mechanism(double r1_, double r2_, double r3_, double r4_, double r5_, double alpha_)
        : r1(r1_), r2(r2_), r3(r3_), r4(r4_), r5(r5_), alpha(alpha_) {
        k1 = r1 / r2;
        k2 = r1 / r4;
        k3 = (r1 * r1 + r2 * r2 + r4 * r4 - r3 * r3) / (2 * r2 * r4);
    }

    double computeA() {
        return (k3 - k2 * cos(theta2)) + (cos(theta2) - k1);
    }

    double computeB() {
        return -2 * sin(theta2);
    }

    double computeC() {
        return (k3 - k2 * cos(theta2)) - (cos(theta2) - k1);
    }

    double computeTheta4() {
        double A = computeA();
        double B = computeB();
        double C = computeC();
        return 2 * atan((-B - sqrt(B * B - 4 * A * C)) / (2 * A));
    }

    double computeTheta3() {
        return asin((r4 * sin(theta4) - r2 * sin(theta2)) / r3);
    }

    void simulate() {
        FILE* fileA = fopen("data/trajectory_A.csv", "w");
        FILE* fileB = fopen("data/trajectory_B.csv", "w");
        FILE* velA = fopen("data/velocity_A.csv", "w");
        FILE* velB = fopen("data/velocity_B.csv", "w");
        FILE* accA = fopen("data/acceleration_A.csv", "w");
        FILE* accB = fopen("data/acceleration_B.csv", "w");

        if (!fileA || !fileB || !velA || !velB || !accA || !accB) {
            std::cerr << "Error opening output files.\n";
            return;
        }

        fprintf(fileA, "xA,yA\n");
        fprintf(fileB, "xB,yB\n");
        fprintf(velA, "theta2,vx,vy\n");
        fprintf(velB, "theta2,vx,vy\n");
        fprintf(accA, "theta2,ax,ay\n");
        fprintf(accB, "theta2,ax,ay\n");

        double xA_vals[N], yA_vals[N], xB_vals[N], yB_vals[N], theta2_vals[N];
        int index = 0;

        for (theta2 = 0; theta2 <= 4 * M_PI && index < N; theta2 += D_THETA2) {
            theta4 = computeTheta4();
            theta3 = computeTheta3();

            double xA = r2 * cos(theta2) + r5 * cos(theta3 + alpha);
            double yA = r2 * sin(theta2) + r5 * sin(theta3 + alpha);
            double xB = r1 + r4 * cos(theta4);
            double yB = r4 * sin(theta4);

            fprintf(fileA, "%lf,%lf\n", xA, yA);
            fprintf(fileB, "%lf,%lf\n", xB, yB);

            xA_vals[index] = xA;
            yA_vals[index] = yA;
            xB_vals[index] = xB;
            yB_vals[index] = yB;
            theta2_vals[index] = theta2;
            ++index;
        }

        const int count = index;

        // Central-difference velocity: v[i] = (x[i+1] - x[i-1]) / (2*D_THETA2)
        for (int i = 1; i < count - 1; ++i) {
            double vxA = (xA_vals[i + 1] - xA_vals[i - 1]) / (2 * D_THETA2);
            double vyA = (yA_vals[i + 1] - yA_vals[i - 1]) / (2 * D_THETA2);
            double vxB = (xB_vals[i + 1] - xB_vals[i - 1]) / (2 * D_THETA2);
            double vyB = (yB_vals[i + 1] - yB_vals[i - 1]) / (2 * D_THETA2);

            fprintf(velA, "%lf,%lf,%lf\n", theta2_vals[i], vxA, vyA);
            fprintf(velB, "%lf,%lf,%lf\n", theta2_vals[i], vxB, vyB);
        }

        // Central-difference acceleration: a[i] = (x[i+1] - 2*x[i] + x[i-1]) / D_THETA2^2
        for (int i = 1; i < count - 1; ++i) {
            double axA = (xA_vals[i + 1] - 2 * xA_vals[i] + xA_vals[i - 1]) / (D_THETA2 * D_THETA2);
            double ayA = (yA_vals[i + 1] - 2 * yA_vals[i] + yA_vals[i - 1]) / (D_THETA2 * D_THETA2);
            double axB = (xB_vals[i + 1] - 2 * xB_vals[i] + xB_vals[i - 1]) / (D_THETA2 * D_THETA2);
            double ayB = (yB_vals[i + 1] - 2 * yB_vals[i] + yB_vals[i - 1]) / (D_THETA2 * D_THETA2);

            fprintf(accA, "%lf,%lf,%lf\n", theta2_vals[i], axA, ayA);
            fprintf(accB, "%lf,%lf,%lf\n", theta2_vals[i], axB, ayB);
        }

        fclose(fileA);
        fclose(fileB);
        fclose(velA);
        fclose(velB);
        fclose(accA);
        fclose(accB);
    }
};

int main() {
    Mechanism mechanism(280, 140, 240, 190, 150, 0.6);
    mechanism.simulate();
    return 0;
}
