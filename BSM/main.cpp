#include <iostream>
#include <cmath>    /* std::erf, std::exp, std::log, std::sqrt */
#include <string>

double norm_cdf(double x);
double bs_call(double S, double K, double σ, double r, double T);
double bs_put(double S, double K, double σ, double r, double T);

int main(int argc, char* argv[]) {
	if (argc != 12) {
		std::cerr << "Usage: " << argv[0]
			<< " -t <C|P> -s <S0> -k <K> -v <σ> -r <r> -T <T>\n";
		return 1;
	}

	char   type = argv[1][0];
	double S    = 0, K = 0, σ = 0, r = 0, T = 0;
	for (int i = 2; i < argc; i += 2) {
		std::string flag = argv[i];
		std::string val  = argv[i+1];
		if      (flag == "-t") type = val[0];
		else if (flag == "-s") S    = std::stod(val);
		else if (flag == "-k") K    = std::stod(val);
		else if (flag == "-v") σ    = std::stod(val);
		else if (flag == "-r") r    = std::stod(val);
		else if (flag == "-T") T    = std::stod(val);
		else {
			std::cerr << "Unknown flag: " << flag << "\n";
			return 1;
		}
	}

	if ((type!='C' && type!='c' && type!='P' && type!='p')
			|| S <= 0 || K <= 0 || σ <= 0 || T <= 0) {
		std::cerr << "Invalid inputs; check Usage above.\n";
		return 1;
	}

	double price = (type=='C' || type=='c')
		? bs_call(S, K, σ, r, T)
		: bs_put (S, K, σ, r, T);

	std::cout << "Price: " << price << "\n";
	return 0;
}

/* cumulative normal distribution via std::erf */
double norm_cdf(double x) {
	return 0.5 * (1.0 + std::erf(x / std::sqrt(2.0)));
}

double bs_call(double S, double K, double σ, double r, double T) {
	double d1 = (std::log(S/K) + (r + 0.5*σ*σ)*T) / (σ*std::sqrt(T));
	double d2 = d1 - σ*std::sqrt(T);
	return S*norm_cdf(d1) - K*std::exp(-r*T)*norm_cdf(d2);
}

double bs_put(double S, double K, double σ, double r, double T) {
	double d1 = (std::log(S/K) + (r + 0.5*σ*σ)*T) / (σ*std::sqrt(T));
	double d2 = d1 - σ*std::sqrt(T);
	return K*std::exp(-r*T)*norm_cdf(-d2) - S*norm_cdf(-d1);
}
