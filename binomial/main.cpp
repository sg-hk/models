#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

int main(int argc, char* argv[]) {
	if (argc != 14) {
		std::cerr << "Usage: " << argv[0]
			<< " <C|P> -n <steps> -s <S0> -k <K> "
			"-r <r> -T <T> -v <sigma>\n";
		return 1;
	}

	char type = argv[1][0];
	int  n    = 0;
	double S0=0, K=0, r=0, T=0, sigma=0;
	for (int i = 2; i < argc; i += 2) {
		std::string flag = argv[i];
		double      val  = std::stod(argv[i+1]);
		if      (flag=="-n") n     = int(val);
		else if (flag=="-s") S0    = val;
		else if (flag=="-k") K     = val;
		else if (flag=="-r") r     = val;
		else if (flag=="-T") T     = val;
		else if (flag=="-v") sigma = val;
		else {
			std::cerr<<"Unknown flag "<<flag<<"\n";
			return 1;
		}
	}
	if (n <= 0 || S0<=0||K<=0||T<=0||sigma<=0) {
		std::cerr<<"All numeric args must be positive\n";
		return 1;
	}

	double dt    = T / n;
	double sig_d = sigma / std::sqrt(n); // per-step vol
	double u     = std::exp(sig_d);
	double d     = std::exp(-sig_d);
	double p     = (std::exp(r*dt) - d) / (u - d);

	// initialize payoff at maturity
	std::vector<double> payoff(n+1);
	for (int j = 0; j <= n; ++j) {
		double ST = S0 * std::pow(u, j) * std::pow(d, n-j);
		payoff[j] = (type=='C' || type=='c')
			? std::max(ST - K, 0.0)
			: std::max(K - ST, 0.0);
	}

	// backward induction
	for (int step = n-1; step >= 0; --step) {
		for (int j = 0; j <= step; ++j) {
			payoff[j] = std::exp(-r*dt) *
				(p*payoff[j+1] + (1-p)*payoff[j]);
		}
	}

	std::cout << "Option price: " << payoff[0] << "\n";
	return 0;
}
