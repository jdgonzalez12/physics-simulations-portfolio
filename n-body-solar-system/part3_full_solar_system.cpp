// part3_full_solar_system.cpp
//
// Generic N-body gravitational simulation, instantiated with the Sun plus
// the 8 planets on circular, coplanar initial orbits. Port of pte3.cpp --
// Proyecto 4, Part 3.
//
// Each body feels the pairwise Newtonian gravitational pull of every other
// body (eq. 6 in the assignment):
//     F_i = sum_{j != i} -G * m_i * m_j / r_ij^2 * r_ij_hat
//
// Same true Velocity Verlet integrator as Parts 1-2 (see
// part1_two_body.cpp for the derivation):
//     a_{n+1} = F / m_i
//     r_{n+1} = r_n + v_n*dt + 0.5*a_n*dt^2
//     v_{n+1} = v_n + 0.5*(a_{n+1} + a_n)*dt
// The velocity update averages the OLD and NEW (recomputed) accelerations
// rather than using only the old one (which would be symplectic Euler).
//
// Simulated duration: the original code had a comment claiming "100 anos"
// (100 years) but actually computed 365.25*50*24*3600 seconds = 50 years,
// which is enough for Jupiter (~11.9-year period) to complete several
// orbits and Saturn (~29.5-year period) a good fraction of one, but not
// enough for Uranus (~84-year period) or Neptune (~165-year period) to
// show more than a small arc. This version simulates the full ~165-year
// Neptune period instead, with dt = 2 days (~44 samples per Mercury orbit,
// ~30,000 rows total, a CSV of a few MB) so every planet, including
// Neptune, completes at least one full orbit while inner-planet orbits
// still look smooth.

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#define PI 3.14159265358979323846
#define G 6.67e-11

class Vector {
public:
    double x, y;

    Vector(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    Vector operator+(const Vector& o) const { return {x + o.x, y + o.y}; }
    Vector operator-(const Vector& o) const { return {x - o.x, y - o.y}; }
    Vector operator*(double s) const { return {x * s, y * s}; }
    Vector operator/(double s) const { return {x / s, y / s}; }
    double norm() const { return std::sqrt(x * x + y * y); }
};

class Body {
public:
    double mass, radius;
    Vector pos, vel, accel;

    Body(double mass_, double radius_, Vector pos_, Vector vel_)
        : mass(mass_), radius(radius_), pos(pos_), vel(vel_), accel(0, 0) {}
};

class SolarSystem {
public:
    std::vector<Body> bodies;

    void add(const Body& b) {
        bodies.push_back(b);
    }

    void computeAccelerations() {
        for (auto& b : bodies)
            b.accel = {0, 0};

        for (size_t i = 0; i < bodies.size(); ++i) {
            for (size_t j = 0; j < bodies.size(); ++j) {
                if (i == j) continue;
                Vector r = bodies[j].pos - bodies[i].pos;
                double d = r.norm();
                bodies[i].accel = bodies[i].accel + r * (G * bodies[j].mass / std::pow(d, 3));
            }
        }
    }

    // Velocity Verlet step for every body. Assumes each body's accel field
    // already holds the acceleration evaluated at the CURRENT position
    // (a_n) -- callers must call computeAccelerations() once before the
    // first step().
    void step(double dt) {
        // Save the "old" accelerations a_n before computeAccelerations()
        // overwrites them below.
        std::vector<Vector> a_old(bodies.size());
        for (size_t i = 0; i < bodies.size(); ++i)
            a_old[i] = bodies[i].accel;

        // (1) Update positions using the current (old) accelerations.
        for (size_t i = 0; i < bodies.size(); ++i) {
            bodies[i].pos = bodies[i].pos + bodies[i].vel * dt + a_old[i] * (0.5 * dt * dt);
        }

        // (2) Recompute accelerations at the new positions -> a_{n+1}.
        computeAccelerations();

        // (3) Update velocities using the AVERAGE of the old and new
        // accelerations (true Velocity Verlet, eq. 4 in the assignment).
        for (size_t i = 0; i < bodies.size(); ++i) {
            bodies[i].vel = bodies[i].vel + (a_old[i] + bodies[i].accel) * (0.5 * dt);
        }
    }

    void simulate(double t_final, double dt, const std::string& filename) {
        std::ofstream file(filename);
        file << "t";
        for (size_t i = 0; i < bodies.size(); ++i)
            file << ",x" << (i + 1) << ",y" << (i + 1);
        file << "\n";

        computeAccelerations(); // initial accelerations (a_0)

        for (double t = 0; t <= t_final; t += dt) {
            file << t;
            for (const auto& b : bodies)
                file << "," << b.pos.x << "," << b.pos.y;
            file << "\n";

            step(dt);
        }

        file.close();
    }
};

int main() {
    double year = 365.25 * 24 * 3600;
    double t_final = 165 * year; // full Neptune orbital period (see header comment)
    double dt = 2 * 24 * 3600;   // 2 days

    SolarSystem system;

    system.add(Body(1.989e30, 6.96e8, {0, 0}, {0, 0}));            // Sun
    system.add(Body(3.3e23, 2.44e6, {5.8e10, 0}, {0, 47870}));     // Mercury
    system.add(Body(4.87e24, 6e6, {1.08e11, 0}, {0, 35000}));      // Venus
    system.add(Body(5.972e24, 6.3e6, {1.5e11, 0}, {0, 29780}));    // Earth
    system.add(Body(6.42e23, 3.39e6, {2.28e11, 0}, {0, 24000}));   // Mars
    system.add(Body(1.9e27, 7e7, {7.78e11, 0}, {0, 13000}));       // Jupiter
    system.add(Body(5.68e26, 5.82e7, {1.43e12, 0}, {0, 9640}));    // Saturn
    system.add(Body(8.68e25, 2.54e7, {2.87e12, 0}, {0, 6800}));    // Uranus
    system.add(Body(1.02e26, 2.4e7, {4.5e12, 0}, {0, 5430}));      // Neptune

    system.simulate(t_final, dt, "data/solar_system.csv");

    std::cout << "Solar system simulation completed.\n";
    return 0;
}
