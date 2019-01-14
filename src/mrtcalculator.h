#ifndef CALCULATEMRT_MRTCALCULATOR_H_
#define CALCULATEMRT_MRTCALCULATOR_H_

#include <vector>

#include <SC_SurfaceCalc.hpp>

#include "date.h"
#include "struct.h"
#include "jsonstruct.h"
#include "fisheye/fisheye.h"
#include "vfcalculator.h"
#include "fractioncalculator.h"

#define M_PI 3.14159265358979323846

namespace calculate {
class MRTCalculator {
public:

	double getSunHours(Date *date, const unsigned char *data, double lat, double lng);
	void sunPositionToImageCoordinates(double dZenithAngle, double dAzimuth, int &i, int &j);
	void getSunPos(const Date &date, double lat, double lng, double &dZenithAngle, double &dAzimuth);
	bool getSunVisible(double dZenithAngle, double dAzimuth, const unsigned char *data);

	void calculate(Date &date, const std::vector<unsigned char *> &data, const file::GeoStruct &geo, const std::vector<double> &viewFactor, std::string *output, const std::vector<std::vector<double>> &factors, const file::JsonStruct *json, std::vector<unsigned char*>&, Fisheye *fisheye, VFCalculator &vfCalculator, std::vector<unsigned char*>& seg_data);

	void getSegShadow(unsigned char* segmented, float* shadow, unsigned char* shaded, unsigned char* non_shaded, vector<int>& time);
	//buffered date to reduce redundant calculation
	Date *date = nullptr;

	std::vector<unsigned char*> shaded_buffer;
	std::vector<unsigned char*> non_shaded_buffer;

	std::vector<unsigned char*> shaded_fe_buffer;
	std::vector<unsigned char*> non_shaded_fe_buffer;

  std::vector<unsigned char*> surf_temp, surf_emission;
  std::vector<unsigned char*> surf_temp_fisheye, surf_emission_fisheye;

	SC_SurfaceCalc *scSuf = nullptr;
  FractionCalculator *fc = nullptr;

	// view factor for shaded and nonshaded fisheye
	std::vector<std::vector<double>> shaded_vf, non_shaded_vf;

	std::vector<std::vector<int>> acc_tem_time;

	int pre_day;
	int time_step;
	/**
		@brief Initialize depth library pointer
	*/
	void init(int timeStep);

	~MRTCalculator() {
		if (date != nullptr) delete date;
		//TODO: if (scSuf != nullptr) delete scSuf;
		for (auto p : shaded_buffer) {
			if (p != nullptr) free(p);
		}

		for (auto p : non_shaded_buffer) {
			if (p != nullptr) free(p);
		}

		for (auto p : shaded_fe_buffer) {
			if (p != nullptr) free(p);
		}

		//shaded_buffer.clear();
		for (auto p : non_shaded_fe_buffer) {
			if (p != nullptr) free(p);
		}
	}

private:
	const double pi = M_PI;
	const double twopi = 2.0 * pi;
	const double rad = pi / 180.0;
	const double dEarthMeanRadius = 6371.01;
	const double dAstronomicalUnit = 149597890;
};
} // calculate

#endif // !CALCULATEMRT_MRTCALCULATOR_H_