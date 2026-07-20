#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>

//                          APPLIED PROGRAMMING - PROJECT 1
//                                     PART 1
//                  Parabolic projectile motion WITHOUT air drag

#define G 9.81

class Projectile {
public:
    double x, y, vx, vy, m, theta, vi;

    Projectile(double x, double y, double vx, double vy, double m)
    : x(x), y(y), vx(vx), vy(vy), m(m) {
        vi = sqrt(vx * vx + vy * vy);  // Initial speed magnitude
        theta = atan2(vy, vx);         // Initial launch angle
    }

    // Advance the projectile's state by one time step dt
    void step(double dt) {
        x += vi * cos(theta) * dt;
        y += vi * sin(theta) * dt - 0.5 * G * dt * dt;
        vx = vi * cos(theta);
        vy = vi * sin(theta) - G * dt;

        theta = atan2(vy, vx);   // Recompute the angle
        vi = sqrt(vx * vx + vy * vy);
    }
};

// Closed-form (analytic) estimate of the maximum height, h_max = V^2 sin^2(theta) / 2g
double hmax_analytic(const Projectile& p) {
    return (p.vy * p.vy) / (2 * G);
}

// Closed-form (analytic) estimate of the range, R = V^2 sin(2 theta) / g
double range_analytic(const Projectile& p) {
    double t_total = 2 * p.vy / G;
    return p.vx * t_total;
}

// Brute-force numerical estimate of the range: step until the projectile
// crosses y = 0 with resolution dt.
double range_numeric(Projectile p, double dt) {
    while (true) {
        p.step(dt);
        if (p.y < 0) break;
    }
    return p.x;
}

// Brute-force numerical estimate of the maximum height: step until the
// vertical velocity turns negative, with resolution dt.
double hmax_numeric(Projectile p, double dt) {
    while (true) {
        p.step(dt);
        if (p.vy < 0) break;
    }
    return p.y;
}

// Tabulate the projectile's trajectory using 100 discretized points
void writeTrajectoryCsv(Projectile p, const std::string& filename) {
    std::ofstream file(filename);
    file << "x,y\n"; // CSV header

    double t_total = 2 * p.vy / G; // Total time of flight
    double dt = t_total / 100;     // Time step

    for (int i = 0; i <= 100; ++i) {
        file << std::fixed << std::setprecision(6) << p.x << "," << p.y << "\n";
        p.step(dt);
    }

    file.close();
}

// Sweep dt over a log-spaced range and, for each resolution, recompute the
// numeric hmax/range and compare against the analytic values. This is what
// Fig. 3 of the assignment asks for: percent error vs. search resolution dt.
void writeErrorSweepCsv(const Projectile& p0, double dt_min, double dt_max,
                         int n_points, const std::string& filename) {
    std::ofstream file(filename);
    file << "dt,error_hmax_pct,error_range_pct\n";

    double hmax_real = hmax_analytic(p0);
    double range_real = range_analytic(p0);

    double log_min = std::log10(dt_min);
    double log_max = std::log10(dt_max);

    for (int i = 0; i < n_points; ++i) {
        double frac = (n_points == 1) ? 0.0 : static_cast<double>(i) / (n_points - 1);
        double dt = std::pow(10.0, log_min + frac * (log_max - log_min));

        double hmax_num = hmax_numeric(p0, dt);
        double range_num = range_numeric(p0, dt);

        double error_hmax = 100 * fabs(hmax_real - hmax_num) / hmax_real;
        double error_range = 100 * fabs(range_real - range_num) / range_real;

        file << std::scientific << std::setprecision(8) << dt << ","
             << error_hmax << "," << error_range << "\n";
    }

    file.close();
}

int main() {
    Projectile p1(0, 0, 3, 3, 0.1);

    // Write CSV with 100 discretized trajectory points
    writeTrajectoryCsv(p1, "data/trayectoria.csv");

    // Analytic (theoretical) values
    double hmax_real = hmax_analytic(p1);
    double range_real = range_analytic(p1);

    // Numeric values at a representative dt = 0.001
    double dt = 0.001;
    Projectile p_temp(0, 0, 3, 3, 0.1);
    double hmax_num = hmax_numeric(p_temp, dt);
    double range_num = range_numeric(p_temp, dt);

    // Percent errors
    double error_hmax = 100 * fabs(hmax_real - hmax_num) / hmax_real;
    double error_range = 100 * fabs(range_real - range_num) / range_real;

    // Console report
    std::cout << "\nAnalytic vs numeric results:\n";
    std::cout << "Analytic max height:  " << std::fixed << std::setprecision(6) << hmax_real << " m\n";
    std::cout << "Numeric max height:   " << std::fixed << std::setprecision(6) << hmax_num << " m\n";
    std::cout << "Analytic range:       " << std::fixed << std::setprecision(6) << range_real << " m\n";
    std::cout << "Numeric range:        " << std::fixed << std::setprecision(6) << range_num << " m\n";

    std::cout << "\nPercent errors (dt = " << dt << "):\n";
    std::cout << "Max height error: " << std::fixed << std::setprecision(6) << error_hmax << "%\n";
    std::cout << "Range error:      " << std::fixed << std::setprecision(6) << error_range << "%\n";

    // Sweep dt from 1e-6 to 0.5 s (log-spaced) and write the error trend to CSV
    writeErrorSweepCsv(p1, 1e-6, 0.5, 40, "data/error.csv");
    std::cout << "\nError-vs-dt sweep written to data/error.csv\n";

    return 0;
}
