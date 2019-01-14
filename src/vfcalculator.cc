#include "vfcalculator.h"

#include<iostream>
#include<vector>
#include<utility>
#include<string>

#define M_PI 3.14159265358979323846

#define RING_SIZE 32
#define EDGE_LENGTH 512
#define HALF_LENGTH 256
#define OFFSET 0.5
#define RING_RADIUS 8.0

namespace calculate {
	VFCalculator::VFCalculator() {
		t.clear();
		t.resize(RING_SIZE, 0);
		for (int i = 0; i < EDGE_LENGTH; i++) {
			for (int j = 0; j < EDGE_LENGTH; j++) {
				double d = sqrt(
					(i + OFFSET - HALF_LENGTH)*(i + OFFSET - HALF_LENGTH)
					+ (j + OFFSET - HALF_LENGTH)*(j + OFFSET - HALF_LENGTH)
				);
				unsigned ringIndex = (unsigned)floor(d / RING_RADIUS);
				if (ringIndex < RING_SIZE) t[ringIndex]++;
			}
		}
		for (unsigned i = 0; i<t.size(); i++) {
			t_total += t[i];
		}

		svf_max.resize(RING_SIZE, 0);
		int n2 = 2 * RING_SIZE;
		for (int i = 0; i < RING_SIZE; i++) {
			svf_max[i] = sin(M_PI*(2 * i + 1) / n2) / t[i];
		}
	}

	void VFCalculator::calculate_P(const unsigned char *data, std::vector<std::vector<int>> & res) {
		for (int i = 0; i < EDGE_LENGTH; i++) {
			for (int j = 0; j < EDGE_LENGTH; j++) {
				double d = sqrt(
					(i + OFFSET - HALF_LENGTH)*(i + OFFSET - HALF_LENGTH)
					+ (j + OFFSET - HALF_LENGTH)*(j + OFFSET - HALF_LENGTH)
				);
				int ringIndex = (int)floor(d / RING_RADIUS);
				if (ringIndex < RING_SIZE) {
					int pixel = (j*EDGE_LENGTH + i);
					if (((int)data[pixel])>0 && ((int)data[pixel])<7) {
						res[(int)data[pixel] - 1][ringIndex]++;
					}
				}
			}
		}
	}

	double VFCalculator::calculate_SVF(
			const std::vector<double> &svf_max, const std::vector<int> &p) {
		double SVF = 0;
		for (int i = 0; i<RING_SIZE; i++) {
			SVF += svf_max[i] * p[i];
		}
		SVF *= M_PI / (2 * RING_SIZE);
		return SVF;
	}

	void VFCalculator::calculate(unsigned char *data, std::vector<double> &viewFactor) {
		std::vector<std::vector<int>> p(6, std::vector<int>(RING_SIZE, 0));
		calculate_P(data, p);
		for (int i = 0; i < 6; i++) {
			double vf = calculate_SVF(svf_max, p[i]);
			viewFactor[i] = vf;
		}
	}

	VFCalculator::~VFCalculator() {

	}
} // calculate