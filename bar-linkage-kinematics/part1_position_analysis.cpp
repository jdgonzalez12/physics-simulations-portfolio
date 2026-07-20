// part1_position_analysis.cpp
//
// Part 1: Position analysis of a 5-bar oscillating linkage mechanism.
// Solves the closed-loop vector equation r2 + r3 = r1 + r4 for theta4 via the
// Weierstrass (half-angle tangent) substitution, producing a quadratic in
// x = tan(theta4/2). theta3 follows directly from substitution back into the
// loop-closure equations. Points A (on the r5 arm) and B (tip of r4) are then
// located in Cartesian coordinates as the driving angle theta2 sweeps through
// a full revolution.
//
// Reference configuration (from the assignment brief):
//   r1=280, r2=140, r3=240, r4=190, r5=150, alpha=0.6 rad
//
// Ported/translated from the original Spanish-language barras1.cpp. Closed-
// form algebraic solve, no numerical bug in this part.

#define _USE_MATH_DEFINES // needed for M_PI on MinGW/MSVC
#include <iostream>
#include <cstdio>
#include <cmath>

#define TOL 1e-6
#define MAX_ITER 1000

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

    // Coefficients of the quadratic A*x^2 + B*x + C = 0, x = tan(theta4/2)
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
        FILE* tableFile = fopen("data/theta2_theta3_table.csv", "w");

        if (!fileA || !fileB || !tableFile) {
            std::cerr << "Error opening output files.\n";
            return;
        }

        fprintf(fileA, "xA,yA\n");
        fprintf(fileB, "xB,yB\n");
        fprintf(tableFile, "theta2_deg,theta3_deg\n");

        std::cout << "theta2(deg)\ttheta3(deg)\n"; // theta2, theta3 reported in degrees

        for (theta2 = 0; theta2 <= 2 * M_PI; theta2 += (2 * M_PI) / 1000) {
            theta4 = computeTheta4();
            theta3 = computeTheta3();

            double xA = r2 * cos(theta2) + r5 * cos(theta3 + alpha);
            double yA = r2 * sin(theta2) + r5 * sin(theta3 + alpha);
            double xB = r1 + r4 * cos(theta4);
            double yB = r4 * sin(theta4);

            fprintf(fileA, "%lf,%lf\n", xA, yA);
            fprintf(fileB, "%lf,%lf\n", xB, yB);

            double theta2_deg = theta2 * 180 / M_PI;
            double theta3_deg = theta3 * 180 / M_PI;
            fprintf(tableFile, "%lf,%lf\n", theta2_deg, theta3_deg);

            std::cout << theta2_deg << "\t" << theta3_deg << "\n";
        }

        fclose(fileA);
        fclose(fileB);
        fclose(tableFile);
    }
};

int main() {
    Mechanism mechanism(280, 140, 240, 190, 150, 0.6);
    mechanism.simulate();
    return 0;
}
