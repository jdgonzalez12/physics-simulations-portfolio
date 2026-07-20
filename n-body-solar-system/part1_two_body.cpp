// part1_two_body.cpp
//
// Two-body gravitational simulation (Newtonian N-body, N = 2).
// Port of pte1.cpp -- Proyecto 4, Part 1 ("Sistema planetario").
//
// Integrates the two-body problem under Newton's law of gravitation,
//     F = -G * m1 * m2 / r^2 * r_hat,
// using true Velocity Verlet integration (as specified by the assignment):
//     a_{n+1} = F / m_i                                  (eq. 2)
//     r_{n+1} = r_n + v_n*dt + 0.5*a_n*dt^2               (eq. 3)
//     v_{n+1} = v_n + 0.5*(a_{n+1} + a_n)*dt              (eq. 4)
//
// Note the velocity update uses the AVERAGE of the OLD acceleration (a_n)
// and the NEW acceleration (a_{n+1}, evaluated at the updated position) --
// not just a_n. Using only a_n for both the position and velocity update
// (v += a*dt) is semi-implicit/symplectic Euler, not Velocity Verlet; it
// still looks like a reasonable orbit but is less accurate and does not
// conserve energy as well over long integrations.
//
// Runs three mass-ratio scenarios (m2/m1_base = 0.5, 1.0, 2.0) to illustrate
// how the trajectory of body 1 depends on the mass ratio between the bodies.

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#define PI 3.14159265358979323846
#define G 6.67e-11
#define TOL 1e-1

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
        // Save the "old" acceleration a_n before it gets overwritten below.
        Vector a1_old = body1.accel;
        Vector a2_old = body2.accel;

        // (1) Update positions using the current (old) accelerations.
        body1.pos = body1.pos + body1.vel * dt + a1_old * (0.5 * dt * dt);
        body2.pos = body2.pos + body2.vel * dt + a2_old * (0.5 * dt * dt);

        // (2) Recompute accelerations at the new positions -> a_{n+1}.
        computeAcceleration();

        // (3) Update velocities using the AVERAGE of the old and new
        // accelerations -- this is the true Velocity Verlet scheme
        // (eq. 4 in the assignment), not plain symplectic Euler.
        body1.vel = body1.vel + (a1_old + body1.accel) * (0.5 * dt);
        body2.vel = body2.vel + (a2_old + body2.accel) * (0.5 * dt);
    }

    // Adaptive step-size search: halves dt until a full step of size dt
    // and two half-steps of size dt/2 agree to within TOL (step-doubling).
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

    // Total mechanical energy (kinetic + gravitational potential), used as
    // a sanity check on the Velocity Verlet fix: a correct implementation
    // should conserve this to a small relative drift over the run.
    double totalEnergy() const {
        double ke = 0.5 * body1.mass * body1.vel.norm() * body1.vel.norm()
                  + 0.5 * body2.mass * body2.vel.norm() * body2.vel.norm();
        double dist = (body1.pos - body2.pos).norm();
        double pe = -G * body1.mass * body2.mass / dist;
        return ke + pe;
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

void simulate_case(double ratio, const std::string& filename) {
    double year = 730.50 * 24 * 60 * 60; // simulated duration (seconds)
    double base_mass = 5.972e24;         // Earth-like mass
    double m1 = 1.989e30;                // Sun-like mass, fixed
    double m2 = base_mass * ratio;
    double distance = 1.5e11;
    double v_orbital = std::sqrt(G * m1 / distance);

    Vector p1 = {0, 0};
    Vector v1 = {0, 0};
    Vector p2 = {distance, 0};
    Vector v2 = {0, v_orbital};

    Body b1(m1, 6.96e8, p1, v1);
    Body b2(m2, 6.3e6, p2, v2);

    TwoBodySystem system(b1, b2);

    double e_initial = system.totalEnergy();
    system.simulate(year, year / 10000, filename);
    double e_final = system.totalEnergy();
    double drift_pct = std::abs((e_final - e_initial) / e_initial) * 100.0;

    std::cout << filename
              << ": E_initial = " << e_initial << " J"
              << ", E_final = " << e_final << " J"
              << ", relative drift = " << drift_pct << " %\n";
}

int main() {
    simulate_case(0.5, "data/two_body_mass_ratio_low.csv");
    simulate_case(1.0, "data/two_body_mass_ratio_equal.csv");
    simulate_case(2.0, "data/two_body_mass_ratio_high.csv");
    std::cout << "Simulations completed.\n";
    return 0;
}
