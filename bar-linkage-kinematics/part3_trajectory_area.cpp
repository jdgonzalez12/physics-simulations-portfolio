// part3_trajectory_area.cpp
//
// Part 3: Enclosed area of point A's closed trajectory over one period of
// theta2, computed with the shoelace formula applied to the polygon traced
// by the discretized (xA, yA) samples.
//
// Ported/translated from the original Spanish-language barras3.cpp. No
// numerical bug in the original math -- this version additionally writes the
// result to data/enclosed_area.csv alongside the existing console print, so
// the value is reproducible without re-running and re-parsing stdout.

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

    // Shoelace-formula area enclosed by the trajectory of point A over one
    // period of theta2, in [0, 2*pi].
    double computeTrajectoryAreaA() {
        double area = 0.0;
        double x1 = 0, y1 = 0;
        bool firstIteration = true;

        for (theta2 = 0; theta2 <= 2 * M_PI; theta2 += (2 * M_PI) / 1000) {
            theta4 = computeTheta4();
            theta3 = computeTheta3();

            double xA = r2 * cos(theta2) + r5 * cos(theta3 + alpha);
            double yA = r2 * sin(theta2) + r5 * sin(theta3 + alpha);

            if (firstIteration) {
                x1 = xA;
                y1 = yA;
                firstIteration = false;
                continue;
            }

            double x2 = xA;
            double y2 = yA;

            area += (x2 - x1) * (y2 + y1) / 2.0;

            x1 = x2;
            y1 = y2;
        }

        return fabs(area);
    }
};

int main() {
    Mechanism mechanism(280, 140, 240, 190, 150, 0.6);
    double area = mechanism.computeTrajectoryAreaA();

    std::cout << "Approximate area enclosed by the trajectory of point A: "
              << area << " square units\n";

    FILE* areaFile = fopen("data/enclosed_area.csv", "w");
    if (!areaFile) {
        std::cerr << "Error opening output file.\n";
        return 1;
    }
    fprintf(areaFile, "area_sq_units,%lf\n", area);
    fclose(areaFile);

    return 0;
}
