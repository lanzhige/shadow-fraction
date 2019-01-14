#ifndef CALCULATEMRT_VFCALCULATOR_H_
#define CALCULATEMRT_VFCALCULATOR_H_

#include<vector>

namespace calculate {
	class VFCalculator {
	public:
		VFCalculator();
		void calculate(unsigned char *data, std::vector<double> &viewFactor);
		~VFCalculator();

		std::vector<double> svf_max;
	private:
		double calculate_SVF(const std::vector<double> &svf_max, const std::vector<int> &p);
		void calculate_P(const unsigned char *data, std::vector<std::vector<int>> &res);
		std::vector<int> t;
		double t_total;
	};
}

#endif // !CALCULATEMRT_VFCALCULATOR_H_
