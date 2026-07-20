#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <cmath>
#include <limits>
#include <vector>
#include <iomanip>
#include <string>
#include <algorithm>

/*                    APPLIED PROGRAMMING - PROJECT 1
                                PART 3
        Ballistic projectile attack & defense with linear air drag      */

// Simulation constants
#define G 9.81
#define B 0.05
#define M 5
#define PI M_PI

// Base is a 500 m wide zone centered at x = 5250 m (assignment Fig. 4)
const double BASE_X = 5250.0;

class Projectile {
public:
    double x, y, vx, vy;

    // Initialize state. invertX mirrors vx, used for the defender (P2), who
    // fires from the right side of the field toward the attacker.
    Projectile(double xi, double yi, double vi, double theta, bool invertX = false) {
        double theta_rad = theta * PI / 180.0; // Convert angle to radians
        x = xi;
        y = yi;
        vx = vi * cos(theta_rad);
        vy = vi * sin(theta_rad);
        if (invertX) vx = -vx;
    }

    // Advance the projectile by one time step dt under linear drag
    void step(double dt) {
        double gamma = B / M;
        x += vx * dt - (0.5 * gamma * vx * dt * dt);
        y += vy * dt - ((M * G + B * vy) * dt * dt / (2 * M));
        vx -= gamma * vx * dt;
        vy -= (G + gamma * vy) * dt;
    }

    bool isAirborne() const {
        return y >= 0;
    }
};

// Distance between two projectiles
double distanceBetween(const Projectile& p1, const Projectile& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Approximate the total flight time of a projectile
double computeFlightTime(Projectile p, double dt) {
    double t = 0;
    while (p.isAirborne()) {
        p.step(dt);
        t += dt;
    }
    return t;
}

// Compute the maximum height reached by a projectile
double computeMaxHeight(Projectile p, double dt) {
    double maxHeight = 0;
    while (p.isAirborne()) {
        p.step(dt);
        if (p.vy < 0) break;  // Stop once it starts descending
        maxHeight = p.y;
    }
    return maxHeight;
}

// Result of searching for the best interception for one attacker configuration
struct InterceptionResult {
    double v1 = 0, theta1 = 0;
    double theta2 = 0;
    double delay = 0;
    double distanceFromBase = -1; // -1 sentinel: no interception found
};

// Find the best interception between the attacker (P1) and the defender (P2):
// search over defender launch angles (1 deg resolution) and firing delays,
// preferring the collision point farthest from the attacker's side of the
// base (i.e. the smallest collision x, which stops the attacker earliest).
InterceptionResult findBestInterception(double vi1, double theta1) {
    InterceptionResult best;
    best.v1 = vi1;
    best.theta1 = theta1;

    double bestX = std::numeric_limits<double>::max();
    double bestCollisionDistance = -1;

    // Attacker P1 starts at (0, 0)
    Projectile p1(0, 0, vi1, theta1);

    const double DT_MIN = 0.00001;
    const double DT_MAX = 0.02;

    // Total flight time of P1 (undisturbed)
    double dt1 = DT_MAX;
    double flightTime1 = computeFlightTime(p1, dt1);

    // Defender P2 starts at x = 6500 (1 km behind the base at x = 5250)
    const double vi2 = 300;
    Projectile p2Base(6500, 0, vi2, 0, true);

    double dt2 = DT_MAX;
    double dt = std::min(dt1, dt2) / 10.0;  // Start with a reasonable dt

    const double ANGLE_MIN = 1;
    const double ANGLE_MAX = 89;

    // Maximum height reached by P2 when fired straight up-range at 89 deg
    Projectile p2MaxHeight(6500, 0, vi2, 89, true);
    double p2MaxAltitude = computeMaxHeight(p2MaxHeight, dt);

    const double CRITICAL_VELOCITY = 600.0;

    // Scale the maximum delay search window with P1's speed
    double maxDelay = (vi1 > CRITICAL_VELOCITY) ? flightTime1 * 3.0 : flightTime1 * 0.5;

    const double COLLISION_DISTANCE = 1.0;
    const double TOLERANCE_DISTANCE = 10.0;

    for (int theta2 = ANGLE_MIN; theta2 <= ANGLE_MAX; theta2++) {
        p2Base = Projectile(6500, 0, vi2, theta2, true);

        for (double delay = 0; delay <= maxDelay; delay += dt * 2) {
            Projectile p1Now = p1;
            Projectile p2Now = p2Base;
            double t = 0;

            // If P1 hasn't reached P2's max altitude yet, push the delay out further
            if (vi1 > CRITICAL_VELOCITY && p1Now.y < p2MaxAltitude) {
                delay += dt * 10;
            }

            // Advance P1 up to the firing delay
            while (t < delay && p1Now.isAirborne()) {
                p1Now.step(dt);
                t += dt;
            }

            // Simulate until collision or until either projectile lands
            while (p1Now.isAirborne() && p2Now.isAirborne()) {
                p1Now.step(dt);
                p2Now.step(dt);
                t += dt;

                double d = distanceBetween(p1Now, p2Now);

                // Refine dt when the projectiles are within 10 m of each other
                if (d < TOLERANCE_DISTANCE) {
                    dt = DT_MIN;
                } else {
                    dt = std::min(DT_MAX, dt * 1.1);
                }

                // Check for a collision
                if (d < COLLISION_DISTANCE) {
                    if (p1Now.x < bestX) {
                        bestX = p1Now.x;
                        bestCollisionDistance = d;
                        best.theta2 = theta2;
                        best.delay = delay;
                        best.distanceFromBase = fabs(p1Now.x - BASE_X);
                    }
                    break;
                }
            }
        }
    }

    if (bestCollisionDistance < 0) {
        best.distanceFromBase = -1; // No interception found for this configuration
    }

    return best;
}

// Print a readable console table and write the results to CSV
void reportResults(const std::vector<InterceptionResult>& results, const std::string& filename) {
    std::ofstream file(filename);
    file << "V1,theta1_deg,theta2_deg,delay_s,distance_from_base_m\n";

    std::cout << "\n-------------------------------------------------------------------\n";
    std::cout << std::left
              << std::setw(10) << "V1 (m/s)"
              << std::setw(14) << "theta1 (deg)"
              << std::setw(14) << "theta2 (deg)"
              << std::setw(12) << "delay (s)"
              << std::setw(18) << "dist. from base (m)" << "\n";
    std::cout << "-------------------------------------------------------------------\n";

    for (const auto& r : results) {
        file << std::fixed << std::setprecision(4)
             << r.v1 << "," << r.theta1 << "," << r.theta2 << ","
             << r.delay << "," << r.distanceFromBase << "\n";

        std::cout << std::left << std::fixed << std::setprecision(2)
                  << std::setw(10) << r.v1
                  << std::setw(14) << r.theta1
                  << std::setw(14) << r.theta2
                  << std::setw(12) << r.delay;
        if (r.distanceFromBase >= 0)
            std::cout << std::setw(18) << r.distanceFromBase << "\n";
        else
            std::cout << std::setw(18) << "NO INTERCEPTION" << "\n";
    }
    std::cout << "-------------------------------------------------------------------\n";

    file.close();
}

int main() {
    struct AttackerCase { double v1, theta1; };

    // The 10 attacker configurations from the assignment (Fig. 5)
    const std::vector<AttackerCase> cases = {
        {324, 68}, {300, 63}, {400, 10}, {480, 7}, {500, 80},
        {250, 45}, {1200, 87}, {270, 30}, {890, 2}, {3000, 89}
    };

    std::vector<InterceptionResult> results;
    results.reserve(cases.size());

    for (const auto& c : cases) {
        std::cout << "\n=== Attacker case: V1 = " << c.v1 << " m/s, theta1 = " << c.theta1 << " deg ===\n";
        results.push_back(findBestInterception(c.v1, c.theta1));
    }

    reportResults(results, "data/interception_results.csv");

    return 0;
}
