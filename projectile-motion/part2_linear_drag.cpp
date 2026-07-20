#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

//                          APPLIED PROGRAMMING - PROJECT 1
//                                     PART 2
//                Parabolic projectile motion WITH linear air drag

#define G 9.81
#define b 0.05 // Air drag coefficient (kg/s), global for this simulation

class Projectile {
private:
    double x, y, vi, theta, m, vx, vy;

public:
    Projectile(double x0, double y0, double vi0, double theta0, double m0)
        : x(x0), y(y0), vi(vi0), theta(theta0), m(m0) {
        // Convert the launch angle to radians for the calculations
        double theta_rad = theta * M_PI / 180.0;
        theta = theta_rad; // Store theta in radians
        // Compute vx and vy from vi and theta
        vx = vi * cos(theta);
        vy = vi * sin(theta);
    }

    // Semi-analytic step using the linear-drag equations (assignment Eqs. 6-9)
    void step(double dt) {
        x += vx * dt - 0.5 * vx * b * dt * dt / m;
        y += vy * dt - 0.5 * dt * dt * (vy * b + m * G) / m;
        vx -= vx * b * dt / m;
        vy -= (G + vy * b / m) * dt;
    }

    bool isAirborne() const {
        return y >= 0;
    }

    // Accessors for the private state
    double getX() const { return x; }
    double getY() const { return y; }
};

// Compute the total time of flight
double computeFlightTime(Projectile p) {
    double flightTime = 0.0;
    double dt = 0.001; // Small time step for high precision

    while (p.isAirborne()) {
        p.step(dt);
        flightTime += dt;
    }

    return flightTime;
}

// Generate the trajectory with exactly 100 points
void writeTrajectoryCsv(Projectile p, const std::string& filename) {
    std::ofstream file(filename);
    file << "X,Y\n"; // CSV header

    double flightTime = computeFlightTime(p);
    double dt = flightTime / 99; // 99 steps -> 100 total points

    for (int i = 0; i < 100; ++i) {
        file << std::fixed << std::setprecision(6) << p.getX() << "," << p.getY() << "\n";
        p.step(dt);
    }

    file.close();
    std::cout << "Estimated total flight time: " << flightTime << " seconds\n";
}

int main() {
    Projectile p1(0, 0, 324, 68.1, 5); // x0, y0, initial speed, angle (deg), mass

    writeTrajectoryCsv(p1, "data/trajectory_drag.csv");

    std::cout << "Data saved to data/trajectory_drag.csv\n";

    return 0;
}
