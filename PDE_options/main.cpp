#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

// Thomas algorithm for tridiagonal system:
// a: sub-diagonal (size n), b: diagonal (n+1), c: super-diagonal (n),
// d: RHS/output (n+1)
void tridiagonal_solve(const std::vector<double>& a,
                       const std::vector<double>& b,
                       const std::vector<double>& c,
                       std::vector<double>& d) {
    int n = b.size();
    std::vector<double> cp(n), dp(n);
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double m = b[i] - a[i-1] * cp[i-1];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - a[i-1] * dp[i-1]) / m;
    }
    d[n-1] = dp[n-1];
    for (int i = n-2; i >= 0; --i) {
        d[i] = dp[i] - cp[i] * d[i+1];
    }
}

int main(int argc, char* argv[]) {
    if (argc != 16) {
        std::cerr << "Usage: " << argv[0]
                  << " <C|P> -s <S0> -k <K> -v <sigma> -r <r>"
                     " -T <T> -M <M> -N <N>\n";
        return 1;
    }

    char   type = argv[1][0];
    double S0=0, K=0, sigma=0, r=0, T=0;
    int    M=0, N=0;
    for (int i = 2; i < argc; i += 2) {
        std::string flag = argv[i];
        std::string val  = argv[i+1];
        if      (flag == "-s") S0    = std::stod(val);
        else if (flag == "-k") K     = std::stod(val);
        else if (flag == "-v") sigma = std::stod(val);
        else if (flag == "-r") r     = std::stod(val);
        else if (flag == "-T") T     = std::stod(val);
        else if (flag == "-M") M     = std::stoi(val);
        else if (flag == "-N") N     = std::stoi(val);
        else {
            std::cerr << "Unknown flag: " << flag << "\n";
            return 1;
        }
    }

    if ((type!='C' && type!='c' && type!='P' && type!='p')
     || S0<=0 || K<=0 || sigma<=0 || r<0 || T<=0 || M<3 || N<3) {
        std::cerr << "Invalid inputs; check Usage above.\n";
        return 1;
    }

    // safe upper stock boundary (~4σ tail)
    double Smax = S0 * std::exp((r + 4*sigma)*T);
    double dS   = Smax / M;
    double dt   = T    / N;
    double theta = 0.5;  // Crank–Nicolson weight

    // allocate grids
    std::vector<double> U_old(M+1), U_new(M+1);
    // terminal payoff
    for (int i = 0; i <= M; ++i) {
        double s = i * dS;
        if (type=='C'||type=='c')
            U_old[i] = std::max(s - K, 0.0);
        else
            U_old[i] = std::max(K - s, 0.0);
    }

    // workspace for tridiagonal
    std::vector<double> a(M), b(M+1), c(M), d(M+1);

    // Rannacher: first two full steps as two half-steps BE (theta=1)
    for (int step = 0; step < 2; ++step) {
        double dt2 = dt / 2;
        for (int half = 0; half < 2; ++half) {
            // assemble BE system
            for (int i = 1; i < M; ++i) {
                double i_f = double(i);
                double alpha = 0.5 * sigma*sigma * i_f*i_f;
                double beta  = 0.5 * r * i_f;
                // coefficients for BE: U^{n+1} on LHS
                a[i-1] = - (alpha - beta) * dt2;
                b[i]   = 1 + (2*alpha + r) * dt2;
                c[i]   = - (alpha + beta) * dt2;
                d[i]   = U_old[i];
            }
            // boundary nodes
            b[0] = 1;   d[0] = (type=='C'||type=='c') ? 0 : K*std::exp(-r*(step*2+half+1)*dt2);
            b[M] = 1;   d[M] = (type=='C'||type=='c')
                            ? (Smax - K*std::exp(-r*(step*2+half+1)*dt2))
                            : 0;

            tridiagonal_solve(a, b, c, d);
            std::copy(d.begin(), d.end(), U_old.begin());
        }
    }

    // remaining N-2 steps: Crank–Nicolson
    for (int step = 2; step < N; ++step) {
        // assemble CN system
        for (int i = 1; i < M; ++i) {
            double i_f = double(i);
            double alpha = 0.5 * sigma*sigma * i_f*i_f;
            double beta  = 0.5 * r * i_f;
            // LHS: theta-weighted
            a[i-1] = - theta * (alpha - beta) * dt;
            b[i]   = 1 + theta * (2*alpha + r) * dt;
            c[i]   = - theta * (alpha + beta) * dt;
            // RHS: (1 - theta)-weighted
            d[i] = (1 - theta) * (alpha - beta) * dt * U_old[i-1]
                 + (1 - (1 - theta) * (2*alpha + r) * dt) * U_old[i]
                 + (1 - theta) * (alpha + beta) * dt * U_old[i+1];
        }
        // boundaries
        b[0] = 1;   d[0] = (type=='C'||type=='c') ? 0 : K*std::exp(-r*(step+1)*dt);
        b[M] = 1;   d[M] = (type=='C'||type=='c')
                        ? (Smax - K*std::exp(-r*(step+1)*dt))
                        : 0;

        tridiagonal_solve(a, b, c, d);
        std::copy(d.begin(), d.end(), U_old.begin());
    }

    // linear interpolate to S0
    double x = S0 / dS;
    int    j = std::min(M-1, int(std::floor(x)));
    double w = x - j;
    double price = (1 - w) * U_old[j] + w * U_old[j+1];
    price = std::max(price, 0.0);

    std::cout << "Price: " << price << "\n";
    return 0;
}
