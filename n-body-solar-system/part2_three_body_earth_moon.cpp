// part2_three_body_earth_moon.cpp
//
// Three-body gravitational simulation: Sun-Earth-Moon. Port of pte2.cpp --
// Proyecto 4, Part 2.
//
// Each body feels the pairwise Newtonian gravitational pull of the other
// two (eq. 5/6 in the assignment):
//     F_i = sum_{j != i} -G * m_i * m_j / r_ij^2 * r_ij_hat
//
// Same true Velocity Verlet integrator as Part 1 (see part1_two_body.cpp):
//     a_{n+1} = F / m_i
//     r_{n+1} = r_n + v_n*dt + 0.5*a_n*dt^2
//     v_{n+1} = v_n + 0.5*(a_{n+1} + a_n)*dt
// The velocity update averages the OLD and NEW (recomputed) accelerations
// rather than using only the old one (which would be symplectic Euler).
//
// Free-function style (no system class), matching the structure of the
// original file. Uses the real Earth-Moon mass ratio (81.3) and realistic
// Sun-Earth / Earth-Moon distances.

#include <iostream>
#include <fstream>
#include <cmath>

#define G 6.67e-11
#define PI 3.14159265358979323846
#define TOL 1e-3

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
    Body(double mass_ = 0, double radius_ = 0, Vector pos_ = {}, Vector vel_ = {})
        : mass(mass_), radius(radius_), pos(pos_), vel(vel_), accel({0, 0}) {}
};

void computeAccelerations(Body& b1, Body& b2, Body& b3) {
    Body* bodies[] = {&b1, &b2, &b3};
    for (int i = 0; i < 3; ++i) {
        bodies[i]->accel = {0, 0};
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            Vector r = bodies[j]->pos - bodies[i]->pos;
            double d = r.norm();
            bodies[i]->accel = bodies[i]->accel + r * (G * bodies[j]->mass / pow(d, 3));
        }
    }
}

// Velocity Verlet step for all three bodies. Assumes each body's accel
// field already holds the acceleration evaluated at the CURRENT position
// (a_n) -- callers must call computeAccelerations() once before the first
// step().
void step(Body& b1, Body& b2, Body& b3, double dt) {
    Body* bodies[] = {&b1, &b2, &b3};

    // Save the "old" accelerations a_n before computeAccelerations()
    // overwrites them below.
    Vector a_old[3];
    for (int i = 0; i < 3; ++i) a_old[i] = bodies[i]->accel;

    // (1) Update positions using the current (old) accelerations.
    for (int i = 0; i < 3; ++i) {
        bodies[i]->pos = bodies[i]->pos + bodies[i]->vel * dt + a_old[i] * (0.5 * dt * dt);
    }

    // (2) Recompute accelerations at the new positions -> a_{n+1}.
    computeAccelerations(b1, b2, b3);

    // (3) Update velocities using the AVERAGE of the old and new
    // accelerations (true Velocity Verlet, eq. 4 in the assignment).
    for (int i = 0; i < 3; ++i) {
        bodies[i]->vel = bodies[i]->vel + (a_old[i] + bodies[i]->accel) * (0.5 * dt);
    }
}

int main() {
    double year = 730.50 * 24 * 3600;
    double dt = year / 10000;
    std::ofstream file("data/sun_earth_moon.csv");
    file << "t,x1,y1,x2,y2,x3,y3\n";

    double sun_mass = 1.989e30;
    double earth_mass = 5.972e24;
    double moon_mass = earth_mass / 81.3;

    double distance = 1.5e11; // Sun-Earth distance
    double earth_orbital_v = std::sqrt(G * sun_mass / distance);

    double moon_distance = 3.84e8; // Earth-Moon distance
    double moon_orbital_v = std::sqrt(G * earth_mass / moon_distance);

    Vector sun_pos = {0, 0};
    Vector sun_vel = {0, 0};

    Vector earth_pos = {distance, 0};
    Vector earth_vel = {0, earth_orbital_v};

    Vector moon_pos = {distance + moon_distance, 0};
    Vector moon_vel = {0, earth_orbital_v + moon_orbital_v};

    Body sun(sun_mass, 6.96e8, sun_pos, sun_vel);
    Body earth(earth_mass, 6.3e6, earth_pos, earth_vel);
    Body moon(moon_mass, 1.74e6, moon_pos, moon_vel);

    computeAccelerations(sun, earth, moon); // initial accelerations (a_0)

    for (double t = 0; t <= year; t += dt) {
        file << t << "," << sun.pos.x << "," << sun.pos.y << ","
             << earth.pos.x << "," << earth.pos.y << ","
             << moon.pos.x << "," << moon.pos.y << "\n";

        step(sun, earth, moon, dt);
    }

    file.close();
    std::cout << "Sun-Earth-Moon simulation completed.\n";
    return 0;
}
