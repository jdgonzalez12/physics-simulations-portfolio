// part1_sun_earth_exact_ratio.cpp
//
// Two-body gravitational simulation, Part 1 special case: exact Sun-Earth
// mass ratio (m_sun / m_earth = 332946). Port of parte1_4.cpp.
//
// Same physics and integrator as part1_two_body.cpp -- see that file for a
// detailed derivation. In short, this uses true Velocity Verlet integration:
//     a_{n+1} = F / m_i
//     r_{n+1} = r_n + v_n*dt + 0.5*a_n*dt^2
//     v_{n+1} = v_n + 0.5*(a_{n+1} + a_n)*dt
// i.e. the velocity update averages the OLD and NEW accelerations, rather
// than using only the old acceleration (which would be symplectic Euler).

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#define PI 3.14159265358979323846
#define G 6.67e-11
#define TOL 1e-4

class Vector {
public:
    double x, y;

    Vector(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}

    Vector operator+(const Vector& o) const { return Vector(x + o.x, y + o.y); }
    Vector operator-(const Vector& o) const { return Vector(x - o.x, y - o.y); }
    Vector operator*(double s) const { return Vector(x * s, y * s); }
    Vector operator/(double s) const { return Vector(x / s, y / s); }

    double norm() const { return std::sqrt(x * x + y * y); }
};

class Body {
public:
    double mass, radius;
    Vector pos, vel, accel;

    Body(double mass_ = 0, double radius_ = 0, Vector pos_ = {}, Vector vel_ = {})
        : mass(mass_), radius(radius_), pos(pos_), vel(vel_), accel({0, 0}) {}
};

class TwoBodySystem {
public:
    Body body1, body2;

    TwoBodySystem(Body b1, Body b2) : body1(b1), body2(b2) {}

    void computeAcceleration() {
        Vector rij = body2.pos - body1.pos;
        double d = rij.norm();
        double d3 = std::pow(d, 3);

        body1.accel = rij * (G * body2.mass / d3);
        body2.accel = rij * (-G * body1.mass / d3);
    }

    // Velocity Verlet step. Assumes body1.accel / body2.accel already hold
    // the acceleration evaluated at the CURRENT position (a_n) -- callers
    // must call computeAcceleration() once before the first step().
    void step(double dt) {
        Vector a1_old = body1.accel;
        Vector a2_old = body2.accel;

        // (1) Update positions using the current (old) accelerations.
        body1.pos = body1.pos + body1.vel * dt + a1_old * (0.5 * dt * dt);
        body2.pos = body2.pos + body2.vel * dt + a2_old * (0.5 * dt * dt);

        // (2) Recompute accelerations at the new positions -> a_{n+1}.
        computeAcceleration();

        // (3) Update velocities using the AVERAGE of the old and new
        // accelerations (true Velocity Verlet, eq. 4 in the assignment).
        body1.vel = body1.vel + (a1_old + body1.accel) * (0.5 * dt);
        body2.vel = body2.vel + (a2_old + body2.accel) * (0.5 * dt);
    }

    // Adaptive step-size search via step-doubling: halves dt until a full
    // step and two half-steps agree to within TOL.
    double findMinTimestep(double dt_initial) {
        Body b1_prev = body1, b2_prev = body2;
        double dt = dt_initial;
        int iter = 0;
        while (true) {
            computeAcceleration();
            step(dt);
            Vector pos_full_b1 = body1.pos;

            body1 = b1_prev;
            body2 = b2_prev;
            computeAcceleration();
            step(dt / 2);
            computeAcceleration();
            step(dt / 2);

            if ((body1.pos - pos_full_b1).norm() < TOL) break;
            if (++iter >= 50) {
                std::cout << "Warning: dt did not converge\n";
                break;
            }
            dt /= 2;
            body1 = b1_prev;
            body2 = b2_prev;
        }
        return dt;
    }

    void simulate(double total_time, double dt_initial, const std::string& output_file) {
        std::ofstream file(output_file);
        file << "t,x1,y1,x2,y2\n";
        for (double t = 0; t <= total_time; t += dt_initial) {
            file << t << "," << body1.pos.x << "," << body1.pos.y << "," << body2.pos.x << "," << body2.pos.y << "\n";
            dt_initial = findMinTimestep(dt_initial);
            computeAcceleration();
            step(dt_initial);
            if ((body1.pos - body2.pos).norm() < body1.radius + body2.radius) {
                std::cout << "Collision detected\n";
                break;
            }
        }
        file.close();
    }
};

int main() {
    double year = 365.25 * 24 * 60 * 60; // 1 year
    double sun_mass = 1.989e30;
    double earth_mass = sun_mass / 332946.0; // exact Sun/Earth mass ratio
    double distance = 1.5e11;
    double v_orbital = std::sqrt(G * sun_mass / distance);

    Vector p1 = {0, 0};
    Vector v1 = {0, 0};
    Vector p2 = {distance, 0};
    Vector v2 = {0, v_orbital};

    Body sun(sun_mass, 6.96e8, p1, v1);
    Body earth(earth_mass, 6.3e6, p2, v2);

    TwoBodySystem system(sun, earth);
    system.simulate(year, year / 1000, "data/sun_earth.csv");

    std::cout << "Sun-Earth simulation completed.\n";
    return 0;
}
