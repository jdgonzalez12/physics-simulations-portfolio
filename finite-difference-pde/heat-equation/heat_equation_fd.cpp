#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

int main() {
    // Problem parameters
    double diffusivity = 110.0; // thermal diffusion coefficient in mm^2/s
    double length = 50.0;       // total rod length in mm
    double total_time = 10.0;   // total simulation time in seconds
    int num_nodes = 20;         // number of spatial nodes

    double dx = length / (num_nodes - 1);
    double dt = 0.5 * dx * dx / diffusivity; // CFL stability criterion

    std::vector<double> temperature(num_nodes, 20.0); // initial temperature vector (20 degrees)
    temperature[0] = 100.0;
    temperature[num_nodes - 1] = 100.0;

    std::vector<double> temperature_prev = temperature;

    std::ofstream file("data/heat_results.csv");
    file << std::fixed << std::setprecision(4);

    // write file header
    file << "time";
    for (int i = 0; i < num_nodes; ++i) {
        file << ",x" << i;
    }
    file << "\n";

    double t_elapsed = 0.0;
    while (t_elapsed < total_time) {
        // save current state to file
        file << t_elapsed;
        for (int i = 0; i < num_nodes; ++i) {
            file << "," << temperature[i];
        }
        file << "\n";

        temperature_prev = temperature; // copy of the data

        for (int i = 1; i < num_nodes - 1; ++i) {
            temperature[i] = dt * diffusivity * (temperature_prev[i - 1] - 2 * temperature_prev[i] + temperature_prev[i + 1]) / (dx * dx) + temperature_prev[i];
        }

        t_elapsed += dt;
        double avg_temperature = 0.0;
        for (double temp : temperature) avg_temperature += temp;
        avg_temperature /= num_nodes;
        std::cout << "t: " << t_elapsed << " s, Average temperature = " << avg_temperature << " C\n";
    }

    file.close();
    std::cout << "Simulation completed and results saved to 'data/heat_results.csv'\n";

    return 0;
}
