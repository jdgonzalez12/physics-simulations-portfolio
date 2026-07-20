#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>

int N = 75;                  // number of spatial points, chosen arbitrarily but reasonably
int M;                       // number of time points

// Computes the maximum stable time step (dt) for the explicit finite-difference
// scheme on the non-uniform spatial grid, following the stability condition
// derived in the accompanying write-up:
//   dtau < h_{i-1} h_i / (r h_{i-1} h_i + sigma^2 x_i^2)
double computeMaxStableDt(const std::vector<double>& S, const std::vector<double>& h, double sigma, double r) {
    double dt_max = 1e10;
    for (size_t i = 1; i < S.size() - 1; i++) { // loop over all spatial points except the endpoints
        double h_i = h[i];
        double h_im1 = h[i - 1];
        double S_i = S[i]; // spatial position
        double denom = r * h_im1 * h_i + sigma * sigma * S_i * S_i;
        double local_dt = (h_im1 * h_i) / denom;
        if (local_dt < dt_max) dt_max = local_dt; // stability condition
    }
    return dt_max;
}

int main() {
    double sigma, r, T, K, Smax;
    std::cout << "Enter the volatility sigma (e.g. 0.2): ";
    std::cin >> sigma;
    std::cout << "Enter the risk-free rate r (e.g. today's rate 0.0436): ";
    std::cin >> r;
    std::cout << "Enter the time to expiration T in years (e.g. 1.0): ";
    std::cin >> T;

    char option_type;
    std::cout << "Option type -- Call or Put? (C/P): ";
    std::cin >> option_type;

    std::cout << "Enter strike price K (e.g. 25): ";
    std::cin >> K;
    std::cout << "Enter max asset price for the grid Smax (e.g. 50): ";
    std::cin >> Smax;

    bool is_call = (option_type == 'C' || option_type == 'c');

    std::vector<double> S(N + 1); // s0,..,sN
    for (int i = 0; i <= N; i++) {
        S[i] = Smax * std::pow(static_cast<double>(i) / N, 2);
    }

    std::vector<double> h(N + 1, 0.0); // filled starting at i = 1, so h[0] is set afterwards
    for (int i = 1; i <= N; i++) {
        h[i] = S[i] - S[i - 1];
        if (h[i] <= 1e-8) h[i] = 1e-8; // avoid overly small spacing that could cause discontinuities
    }
    h[0] = h[1]; // the first value of h, not computed in the loop above (which starts at i = 1), is set equal to the second value.

    double dt = computeMaxStableDt(S, h, sigma, r);
    M = static_cast<int>(T / dt) + 1; // number of time steps M is now known
    dt = T / M;

    std::vector<std::vector<double>> V(M + 1, std::vector<double>(N + 1, 0.0));

    // Terminal payoff at maturity (t = T), branched by option type
    for (int i = 0; i <= N; i++) {
        if (is_call) {
            V[M][i] = std::max(S[i] - K, 0.0); // call payoff: max(S - K, 0)
        } else {
            V[M][i] = std::max(K - S[i], 0.0); // put payoff: max(K - S, 0)
        }
    }

    // Boundary conditions, consistent with the option type above
    for (int j = 0; j <= M; j++) {
        double t = j * dt;
        if (is_call) {
            V[j][0] = 0.0;                                   // call is worthless when S = 0
            V[j][N] = Smax - K * std::exp(-r * (T - t));      // call deep in the money at S = Smax
        } else {
            V[j][0] = K * std::exp(-r * (T - t));             // put value at S = 0, from solving the reduced ODE
            V[j][N] = 0.0;                                    // put is worthless far out of the money at S = Smax
        }
    }


    for (int j = M; j > 0; j--) {
        for (int i = 1; i < N; i++) {
            double h_i = h[i];
            double h_im1 = h[i - 1];
            double S_i = S[i];

            double a = dt * (sigma * sigma * S_i * S_i - r * S_i * h_im1) / (h_im1 * (h_im1 + h_i));
            double b = 1 - r * dt - dt * sigma * sigma * S_i * S_i / (h_im1 * h_i);
            double c = dt * (sigma * sigma * S_i * S_i + r * S_i * h_i) / (h_i * (h_im1 + h_i));

            if (std::isnan(a) || std::isnan(b) || std::isnan(c)) {
                std::cerr << "Error: coefficients a, b, c are NaN at i = " << i << std::endl;
                V[j - 1][i] = 0.0; // cancel the contribution of that NaN grid point
            } else {
                V[j - 1][i] = a * V[j][i - 1] + b * V[j][i] + c * V[j][i + 1]; // linear combination a*V_j,i-1 + b*V_i,j + c*V_i,j+1
            }
        }
    }

    std::string suffix = is_call ? "call" : "put";

    std::ofstream file1("data/grid_V_ij_" + suffix + ".csv");
    for (int j = 0; j <= M; j++) {
        for (int i = 0; i <= N; i++) { // i = 1 to N-1 excludes the endpoints, since those are fixed by the boundary conditions
            file1 << V[j][i];
            if (i < N) file1 << ",";
        }
        file1 << "\n";
    }
    file1.close();

    std::ofstream file2("data/grid_V_S_t_" + suffix + ".csv");
    file2 << "t(days),";
    for (int i = 0; i <= N; i++) {
        file2 << S[i];
        if (i < N) file2 << ",";
    }
    file2 << "\n";
    for (int j = 0; j <= M; j++) {
        double t_days = j * dt * 365; // convert to days for the plotting axis, instead of raw spatial points
        file2 << t_days << ",";
        for (int i = 0; i <= N; i++) {
            file2 << V[j][i];
            if (i < N) file2 << ",";
        }
        file2 << "\n";
    }
    file2.close();



    std::cout << "\nValues of V(S,t=0) across the option's lifetime (S in [0, " << Smax << "]):\n";
    for (int h = 0; h <= 50; h++) {
        double frac = static_cast<double>(h) / 50.0;
        int i = static_cast<int>(frac * N);
        std::cout << "V(S,t=0) at different possible levels of today's asset price, level: " << h << ": S = " << std::fixed << std::setprecision(2) << S[i]
                  << ", V = " << V[0][i] << "\n";
    }

    std::cout << "\nCSV export complete" << std::endl;
    return 0;
}
