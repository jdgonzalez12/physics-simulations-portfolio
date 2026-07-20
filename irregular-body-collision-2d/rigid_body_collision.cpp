// rigid_body_collision.cpp
//
// 2D rigid-body physics engine for an arbitrary irregular shape (loaded from
// an ASCII grid) falling under gravity and bouncing off a flat ground plane,
// exchanging energy between translational and rotational motion on impact.
//
// Physics summary (see README.md for the full derivation):
//   - Mass properties (center of mass, moment of inertia) are computed by
//     discretizing the shape into equal point masses on the input grid.
//   - Between collisions the body moves as a free projectile (Part 1-3
//     baseline) or, optionally, with linear air drag (Part 4).
//   - On ground contact, a closed-form impulsive collision response solves
//     simultaneous conservation of energy and angular momentum about the
//     contact point for the post-collision vertical velocity and angular
//     velocity, with a restitution coefficient e (e=1 elastic for Parts 1-3,
//     e<1 inelastic for Part 4).
//
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include "libbmp.h"

// Gravity, ground height, and drag are expressed in "grid units" (the same
// units as the input shape.txt grid), not SI units -- this is a stylized,
// artistic-pacing choice for a clean visualization, not a claim of physical
// realism. See README.md for the reasoning.
constexpr double GRAVITY = -40.0;     // grid-units / s^2
constexpr double GROUND_Y = 4.0;      // grid-units, height of the ground plane
constexpr double AIR_DRAG_B = 0.05;   // used only by stepWithDrag() (Part 4)

// Rendering only: physics stays entirely in grid-unit space; pixels only
// matter when a frame is rasterized to a BMP file.
constexpr double SCALE_PX = 5.0;      // pixels per grid unit
constexpr int CANVAS_WIDTH = 480;
constexpr int CANVAS_HEIGHT = 270;

struct Vector2D {
    double x = 0.0;
    double y = 0.0;
};

class RigidBody2D {
public:
    RigidBody2D(const std::string& grid_file, double mass_)
        : mass(mass_)
    {
        loadGrid(grid_file);
        computeCenterOfMass();
        computeMomentOfInertia();
        computeContour();
    }

    Vector2D velocity{0.0, 0.0};
    double angular_velocity = 0.0;

    double momentOfInertia() const { return moment_of_inertia; }
    Vector2D centerOfMass() const { return center_of_mass; }

    // Places the body's center of mass at (x, y), keeping its shape/orientation.
    void moveTo(double x, double y) {
        translate(x - center_of_mass.x, y - center_of_mass.y);
    }

    // Part 1-3 baseline: frictionless projectile motion under gravity only.
    void step(double dt) {
        rotate(angular_velocity * dt);
        double dx = velocity.x * dt;
        double dy = velocity.y * dt + 0.5 * GRAVITY * dt * dt;
        translate(dx, dy);
        velocity.y += GRAVITY * dt;
    }

    // Part 4: adds linear air drag (same semi-analytic per-step formula used
    // in projectile-motion/part2_linear_drag.cpp) on top of gravity.
    void stepWithDrag(double dt) {
        rotate(angular_velocity * dt);
        double dx = velocity.x * dt - 0.5 * (AIR_DRAG_B / mass) * velocity.x * dt * dt;
        double dy = velocity.y * dt + 0.5 * (GRAVITY - (velocity.y * AIR_DRAG_B / mass)) * dt * dt;
        translate(dx, dy);
        velocity.y += (GRAVITY - (velocity.y * AIR_DRAG_B / mass)) * dt;
        velocity.x += -(AIR_DRAG_B / mass) * velocity.x * dt;
    }

    // Detects ground contact (any contour point at or below the ground
    // plane) and, if the contact point is still moving into the ground,
    // resolves the collision with a closed-form impulsive response derived
    // from simultaneous energy conservation and angular-impulse/torque
    // balance (Part 2's Eqs. 10-11; see README.md for the derivation).
    // Returns true if a collision was resolved this frame.
    bool resolveGroundCollision(double restitution) {
        int contact_idx = -1;
        double deepest_y = GROUND_Y;
        for (size_t i = 0; i < contour.size(); ++i) {
            if (contour[i].y < deepest_y) {
                deepest_y = contour[i].y;
                contact_idx = static_cast<int>(i);
            }
        }
        if (contact_idx < 0) return false; // no contour point below ground

        // Horizontal moment arm of the contact point about the center of
        // mass -- this is the "r*sin(phi)" term from the PDF's Eq. 11, since
        // only a vertical (normal) impulse is applied here.
        const double rx = contour[contact_idx].x - center_of_mass.x;

        // Vertical velocity of the material point in contact with the
        // ground, accounting for both translation and rotation.
        const double contact_vy = velocity.y + angular_velocity * rx;

        if (contact_vy >= 0.0) return false; // already separating, no impulse needed

        // Frictionless normal impulse J satisfying contact_vy_after = -e * contact_vy_before.
        // Derivation: Δvy = J/mass, Δω = J*rx/I (this is exactly Eq. 11:
        // I*Δω/(r*sinφ) = mass*Δvy). Requiring the contact point's velocity
        // to reverse with restitution e is algebraically equivalent to
        // Eq. 10's energy statement when e = 1 (perfectly elastic).
        const double J = -(1.0 + restitution) * contact_vy /
                          (1.0 / mass + (rx * rx) / moment_of_inertia);

        velocity.y += J / mass;
        angular_velocity += J * rx / moment_of_inertia;

        // Push the body back out of the ground so it doesn't visibly sink in
        // or re-trigger a spurious collision on the very next frame.
        const double penetration = GROUND_Y - contour[contact_idx].y;
        if (penetration > 0.0) translate(0.0, penetration);

        return true;
    }

    double kineticEnergy() const {
        return 0.5 * mass * (velocity.x * velocity.x + velocity.y * velocity.y)
             + 0.5 * moment_of_inertia * angular_velocity * angular_velocity;
    }

    double potentialEnergy() const {
        return -mass * GRAVITY * center_of_mass.y; // GRAVITY < 0, so this is +mass*g*height
    }

    void saveFrameBMP(const std::string& filename) const {
        BmpImg img(CANVAS_WIDTH, CANVAS_HEIGHT);

        for (int y = 0; y < CANVAS_HEIGHT; ++y)
            for (int x = 0; x < CANVAS_WIDTH; ++x)
                img.set_pixel(x, y, 255, 255, 255);

        // Draw the ground line for visual reference.
        int ground_px = CANVAS_HEIGHT - static_cast<int>(std::round(GROUND_Y * SCALE_PX));
        if (ground_px >= 0 && ground_px < CANVAS_HEIGHT) {
            for (int x = 0; x < CANVAS_WIDTH; ++x) img.set_pixel(x, ground_px, 190, 190, 190);
        }

        // Draw each contour point as a filled SCALE_PX x SCALE_PX block (its
        // footprint as a grid cell) rather than a single pixel, so the
        // outline reads as a solid silhouette instead of a sparse dotted line.
        const int block = static_cast<int>(std::round(SCALE_PX));
        for (const auto& p : contour) {
            int x0 = static_cast<int>(std::round(p.x * SCALE_PX));
            int y0 = CANVAS_HEIGHT - static_cast<int>(std::round(p.y * SCALE_PX));
            for (int dy = 0; dy < block; ++dy) {
                for (int dx = 0; dx < block; ++dx) {
                    int x = x0 + dx, y = y0 - dy;
                    if (x >= 0 && x < CANVAS_WIDTH && y >= 0 && y < CANVAS_HEIGHT) {
                        img.set_pixel(x, y, 20, 20, 20);
                    }
                }
            }
        }

        img.write(filename);
    }

private:
    double mass;
    double points_with_mass = 0.0;
    double moment_of_inertia = 0.0;
    Vector2D center_of_mass;
    std::vector<Vector2D> contour;
    std::vector<std::string> grid;
    int rows = 0, cols = 0;

    void loadGrid(const std::string& filename) {
        std::ifstream file(filename);
        std::string line;
        while (std::getline(file, line)) grid.push_back(line);
        rows = static_cast<int>(grid.size());
        cols = grid.empty() ? 0 : static_cast<int>(grid[0].size());
    }

    // Row 0 of the text file is the TOP of the drawn shape: grid row r maps
    // to physics y = (rows - 1 - r), so shape.txt reads naturally top-to-bottom.
    double rowToY(int r) const { return static_cast<double>(rows - 1 - r); }

    void computeCenterOfMass() {
        points_with_mass = 0.0;
        center_of_mass = {0.0, 0.0};
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (grid[r][c] != '0') {
                    center_of_mass.x += c;
                    center_of_mass.y += rowToY(r);
                    points_with_mass += 1.0;
                }
            }
        }
        center_of_mass.x /= points_with_mass;
        center_of_mass.y /= points_with_mass;
    }

    void computeMomentOfInertia() {
        moment_of_inertia = 0.0;
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (grid[r][c] != '0') {
                    double dx = c - center_of_mass.x;
                    double dy = rowToY(r) - center_of_mass.y;
                    moment_of_inertia += dx * dx + dy * dy;
                }
            }
        }
        moment_of_inertia *= mass / points_with_mass;
    }

    // Outline ("shell") points: filled cells with at least one background
    // (or out-of-grid) neighbor in their 8-neighborhood. Treating
    // out-of-range neighbors as background ensures cells touching the edge
    // of the input grid are still correctly detected as contour points.
    void computeContour() {
        contour.clear();
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                if (grid[r][c] == '0') continue;
                bool is_edge = false;
                for (int dr = -1; dr <= 1 && !is_edge; ++dr) {
                    for (int dc = -1; dc <= 1 && !is_edge; ++dc) {
                        int nr = r + dr, nc = c + dc;
                        char neighbor = (nr < 0 || nr >= rows || nc < 0 || nc >= cols)
                                        ? '0' : grid[nr][nc];
                        if (neighbor == '0') is_edge = true;
                    }
                }
                if (is_edge) contour.push_back({static_cast<double>(c), rowToY(r)});
            }
        }
    }

    void rotate(double theta) {
        for (auto& p : contour) { p.x -= center_of_mass.x; p.y -= center_of_mass.y; }
        for (auto& p : contour) {
            double x_new = p.x * std::cos(theta) - p.y * std::sin(theta);
            double y_new = p.x * std::sin(theta) + p.y * std::cos(theta);
            p.x = x_new; p.y = y_new;
        }
        for (auto& p : contour) { p.x += center_of_mass.x; p.y += center_of_mass.y; }
    }

    void translate(double dx, double dy) {
        for (auto& p : contour) { p.x += dx; p.y += dy; }
        center_of_mass.x += dx;
        center_of_mass.y += dy;
    }
};

struct RunConfig {
    std::string label;
    bool use_drag;
    double restitution;
    double duration_s;
};

static void runSimulation(const RunConfig& cfg, int fps) {
    RigidBody2D body("shape.txt", 3.3);
    body.moveTo(20.0, 45.0);
    body.velocity = {6.0, 0.0};
    body.angular_velocity = 0.0;

    const double frame_dt = 1.0 / fps;
    const int total_frames = static_cast<int>(cfg.duration_s * fps);

    // Physics substepping: the ground-plane fall speed is large enough
    // (v = sqrt(2*|GRAVITY|*drop_height) ~ 50-60 grid-units/s) that a single
    // physics step per rendered frame lets the body overshoot noticeably
    // past the ground plane before the collision is detected, injecting
    // spurious energy into the elastic (e=1) case. Running several smaller
    // physics substeps per rendered frame keeps that overshoot -- and the
    // resulting energy drift -- small, while the video still renders at a
    // smooth `fps`.
    const int substeps = 48;
    const double dt = frame_dt / substeps;

    std::ofstream energy_log("data/energy_" + cfg.label + ".csv");
    energy_log << "t,kinetic_energy,potential_energy,total_energy,angular_velocity,com_y\n";

    int collision_count = 0;
    double initial_energy = body.kineticEnergy() + body.potentialEnergy();

    for (int frame = 0; frame <= total_frames; ++frame) {
        std::ostringstream frame_name;
        frame_name << "frames/" << cfg.label << "_frame_"
                   << std::setw(5) << std::setfill('0') << frame << ".bmp";
        body.saveFrameBMP(frame_name.str());

        double t = frame * frame_dt;
        double ke = body.kineticEnergy();
        double pe = body.potentialEnergy();
        energy_log << t << "," << ke << "," << pe << "," << (ke + pe) << ","
                   << body.angular_velocity << "," << body.centerOfMass().y << "\n";

        for (int s = 0; s < substeps; ++s) {
            if (cfg.use_drag) body.stepWithDrag(dt);
            else body.step(dt);

            if (body.resolveGroundCollision(cfg.restitution)) ++collision_count;
        }
    }

    double final_energy = body.kineticEnergy() + body.potentialEnergy();
    double drift_pct = 100.0 * std::fabs(final_energy - initial_energy) / std::fabs(initial_energy);

    std::cout << "[" << cfg.label << "] frames=" << (total_frames + 1)
              << " collisions=" << collision_count
              << " initial_energy=" << initial_energy
              << " final_energy=" << final_energy
              << " drift_pct=" << drift_pct << "%\n";
}

int main() {
    const int fps = 30;

    // Parts 1-3: elastic, frictionless baseline (e = 1.0 conserves mechanical
    // energy across bounces, per Eq. 10; this is the numerical validation
    // check reported in README.md).
    RunConfig elastic{"elastic", false, 1.0, 6.0};

    // Part 4: inelastic collisions (e < 1) plus linear air drag -- bounce
    // height and total energy should visibly decay over time.
    RunConfig inelastic{"inelastic", true, 0.75, 8.0};

    runSimulation(elastic, fps);
    runSimulation(inelastic, fps);

    std::cout << "Simulation complete. Frames written to frames/, energy logs to data/.\n";
    std::cout << "Run ffmpeg (see README.md) to assemble frames into a video.\n";
    return 0;
}
