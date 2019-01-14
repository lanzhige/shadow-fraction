#ifndef CALCULATEMRT_MRTCALCULATOR_H_
#define CALCULATEMRT_MRTCALCULATOR_H_

#include <vector>

#include <SC_SurfaceCalc.hpp>

#include "date.h"
#include "struct.h"
#include "jsonstruct.h"
#include "fisheye/fisheye.h"
#include "vfcalculator.h"

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

	double vectorAngle(double x, double y, double z, double azimuth, int direction);

private:
	void calculateUTCI(double airTemperature, double windVelocity, double mrt, double vaporPressure, std::stringstream &ss);
	void initTotalTemperature(Date &date, double lat, double lng, std::vector<unsigned char *> &seg, std::vector<unsigned char *> &depth);

#ifdef _DEBUG
	unsigned char *debug_image_buffer = nullptr;
	std::string debug_postfix[6] = {
		"_0_90.jpg",   //up
		"_0_270.jpg",  //down
		"_0_0.jpg",    //north
		"_180_0.jpg",  //south
		"_90_0.jpg",   //east
		"_270_0.jpg"  //west
	};
	std::vector<std::vector<unsigned char>> color_map = {
		{ 255, 255, 255, 255 },
		{ 0, 0, 255, 255 },
		{ 0, 255, 0, 255 },
		{ 192, 192, 192, 255 },
		{ 128, 128, 128, 255 },
		{ 255,255,0, 255 },
		{ 255, 0, 255, 255 }
	};

	const std::string debug_temp[6] = {
		"_90_0.jpg",
		"_270_0.jpg",
		"_0_90.jpg",
		"_0_270.jpg",
		"_0_0.jpg",
		"_180_0.jpg"
	};

	/**
	@brief Save images for debugging purpose
	@param data Unsigned char pointer of image input
	@param location Directory to output file
	@param size of image data
	*/
	void saveImage(const unsigned char *data, std::string location, int size);
	/**
	@brief Save accumulated temperature images for debugging purpose
	@param data Vector of accumulated temperature input
	@param location Directory to output file
	*/
	void saveTemperatureImage(const vector<int> &data, std::string location);
#endif //_DEBUG

	const double pi = M_PI;
	const double twopi = 2.0 * pi;
	const double rad = pi / 180.0;
	const double dEarthMeanRadius = 6371.01;
	const double dAstronomicalUnit = 149597890;

	const double C_STEFAN_BOLTZMANN = 5.670367 / 100000000;
	const double C_ENV_EMISSION = 0.95;
	const double C_BODY_EMISSION = 0.97;
	const double C_BODY_ABSORPTION = 0.7;
	const double C_LINKE_TURBIDITY = 4.2;
	const double C_SOLAR = 1367;
	const double C_ALBEDO_ENV = 0.3;
	const double C_BOWEN_RATIO = 1;
	const double C_PRESSURE = 101325;

	const double C_BOWEN_RATIO_SKY = 20;
	const double C_BOWEN_RATIO_TREE = 20;
	const double C_BOWEN_RATIO_BUILDING = 20;
	const double C_BOWEN_RATIO_IMPERVIOUS = 20;
	const double C_BOWEN_RATIO_PERVIOUS = 0.3;
	const double C_BOWEN_RATIO_MOVINGOBJECTS = 20;

	double getVaporPressure(double relativeHumidity, double airTemperature) {
		return relativeHumidity / 100 * (6.11 * pow(10, 7.5 * airTemperature / (237.3 + airTemperature)));
	}

	// Barometric formula (outpt)
	double getAtmosphericPressure(double altitude) {
		return 101325 * pow(1 - 2.25577 / 100000 * altitude, 5.25588);
	}

	// VDI 3789 Appendix B (B1)
	double getRelativeOpticalAirMass(double gamma) {
		if (gamma<0) return 0.0408;
		return 1.0 / (sin(gamma * M_PI / 180) + 0.50572 * pow(gamma + 6.07995, -1.6364));
	}

	// VDI 3789 Appendix B (B2)
	double getOpticalThickness(double relativeOpticalAirMass) {
		return 1.0 / (0.9 * relativeOpticalAirMass + 9.4);
	}

	// VDI 3789 3.1.1 (10, 11, 12, 13)
	double getExtraterrestialRadiation(int dayOfYear) {
		return C_SOLAR * (1 + 0.03344 * cos((dayOfYear * 0.9856 - 2.72) * M_PI / 180));
	}

	double getDirectRadiation(double extraterrestialRadiation, double gamma, double opticalThickness, double relativeOpticalAirMass, double atmosphericPressure, double cloudCover) {
		if (gamma<0) return 0;
		return extraterrestialRadiation * sin(gamma * M_PI / 180)
			* exp(-C_LINKE_TURBIDITY * opticalThickness * relativeOpticalAirMass * atmosphericPressure / C_PRESSURE)
			* (1 - cloudCover / 8);
	}

	double getDirectRadiationPerpendicular(double extraterrestialRadiation, double gamma, double phi, double opticalThickness, double relativeOpticalAirMass, double atmosphericPressure, double cloudCover) {
		return extraterrestialRadiation * getCosMu(phi, 90 - gamma, phi, gamma)
			* exp(-C_LINKE_TURBIDITY * opticalThickness * relativeOpticalAirMass * atmosphericPressure / C_PRESSURE)
			* (1 - cloudCover / 8);
	}

	// VDI 3789 (19, 20) (output)
	double getGlobalRadiation(double extraterrestialRadiation, double gamma, double atmosphericPressure, double cloudCover) {
		if (gamma<0) return 0;
		double sinGamma = sin(gamma * M_PI / 180);
		return 0.84 * extraterrestialRadiation * sinGamma
			* exp(-0.027 * atmosphericPressure / C_PRESSURE * C_LINKE_TURBIDITY / sinGamma)
			* (1 - 0.72 * pow(cloudCover / 8, 3.2));
	}

	double getTransmittanceOfDirectRadiation(double directRadiationCloudless, double extraterrestialRadiation) {
		return directRadiationCloudless / extraterrestialRadiation;
	}

	double getDiffuseRadiationCloudless(double globalRadiationCloudless, double directRadiationCloudless, double transmittanceOfDirectRadiation, double skyViewFactor, double sunVisible) {
		double res = getDiffuseRadiationIsotropic(globalRadiationCloudless, directRadiationCloudless, transmittanceOfDirectRadiation, skyViewFactor);
		if (sunVisible)
			res += getDiffuseRadiationAnisotropic(globalRadiationCloudless, directRadiationCloudless, transmittanceOfDirectRadiation);
		return res;
	}

	// VDI 3789 (24d, 24e, 27) with Omega_H = skyViewFactor
	double getDiffuseRadiationOvercast(double globalRadiationCloudless, double skyViewFactor) {
		return 0.28 * globalRadiationCloudless * skyViewFactor;
	}

	double getDiffuseRadiation(double diffuseRadiationCloudless, double diffuseRadiationOvercast, double cloudCover) {
		return diffuseRadiationCloudless * (1 - cloudCover / 8) + diffuseRadiationOvercast * cloudCover / 8;
	}

	double getDiffuseRadiationIsotropic(double globalRadiationCloudless, double directRadiationCloudless, double transmittanceOfDirectRadiation, double skyViewFactor) {
		return (globalRadiationCloudless - directRadiationCloudless) * (1 - transmittanceOfDirectRadiation) * skyViewFactor;
	}

	// VDI 3789 (21, 22, 28) with beta = 0
	double getDiffuseRadiationAnisotropic(double globalRadiationCloudless, double directRadiationCloudless, double transmittanceOfDirectRadiation) {
		return (globalRadiationCloudless - directRadiationCloudless) * transmittanceOfDirectRadiation;
	}

	double getAtmosphericRadiation(double airTemperature, double vaporPressure, double cloudCover) {
		return C_STEFAN_BOLTZMANN * pow(airTemperature + 273.15, 4) * (0.82 - 0.25 * pow(10, -0.0945 * vaporPressure))
			* (1 + 0.21*pow(cloudCover / 8, 2.5));
	}

	std::pair<double, double> getSurfaceTemperatureAndLongWaveRadiation(double airTemperature, double atmosphericRadiation, double globalRadiation, double windVelocity, const std::string& fractionType) {
		double T_k = airTemperature + 273.15; // air temperature in Kelvin
		double T_s = T_k; // initialize surface temperature with air temperature
		double E, Q, B;

		if (fractionType == "TREE") {
			E = C_ENV_EMISSION * C_STEFAN_BOLTZMANN * pow(T_s, 4) + (1 - C_ENV_EMISSION)*atmosphericRadiation;
			T_s = T_k;
			return{ T_s - 273.15, E };
		}

		for (int i = 0; i < 1; i++) {
			// VDI 3789 (43)
			E = C_ENV_EMISSION * C_STEFAN_BOLTZMANN * pow(T_s, 4) + (1 - C_ENV_EMISSION)*atmosphericRadiation;
			
			Q = C_BODY_ABSORPTION * globalRadiation + atmosphericRadiation - E;
			if (fractionType == "SKY") {
				B = Q > 0 ? -0.19*Q : -0.32*Q;
				T_s = T_k + (Q + B) / ((6.2 + 4.26 * windVelocity) * (1 + 1 / C_BOWEN_RATIO));
			}
			else if (fractionType == "BUILDING") {
				B = -0.51*Q - 34;
				T_s = T_k + (Q + B) / ((6.2 + 4.26 * windVelocity) * (1 + 1 / C_BOWEN_RATIO_BUILDING));
			}
			else if (fractionType == "IMPERVIOUS") {
				B = -0.70*Q - 38;
				T_s = T_k + (Q + B) / ((6.2 + 4.26 * windVelocity) * (1 + 1 / C_BOWEN_RATIO_IMPERVIOUS));
			}
			else if (fractionType == "PERVIOUS") {
				B = -0.34*Q - 31;
				T_s = T_k + (Q + B) / ((6.2 + 4.26 * windVelocity) * (1 + 1 / C_BOWEN_RATIO_PERVIOUS));
			}
			else if (fractionType == "MOVINGOBJECTS") {
				B = -0.34*Q - 31;
				T_s = T_k + (Q + B) / ((6.2 + 4.26 * windVelocity) * (1 + 1 / C_BOWEN_RATIO_PERVIOUS));
				//B = Q > 0 ? -0.19*Q : -0.32*Q;
				//T_s = T_k + (Q + B) / ((6.2 + 4.26 * windVelocity) * (1 + 1 / C_BOWEN_RATIO_MOVINGOBJECTS));
			}
		}
		//T_s = 42+273.15;
		E = C_ENV_EMISSION * C_STEFAN_BOLTZMANN * pow(T_s, 4) + (1 - C_ENV_EMISSION)*atmosphericRadiation;
		return{ T_s - 273.15, E };
		//return [42.0,E];
	}

  void getSurfaceTemperature(double &Ts, double &E, double A, double G, double wind, int type) {
    std::cout << "A: " << A << " G: " << G << std::endl;
    double Tk = Ts;
    double Q, B;
    if (type == 2) {
      E = C_ENV_EMISSION * C_STEFAN_BOLTZMANN * pow(Ts, 4) + (1 - C_ENV_EMISSION) * A;
      Ts = Tk - 273.15;
      return;
    }

    for (int i = 0; i < 1000; i++) {
      E = C_ENV_EMISSION * C_STEFAN_BOLTZMANN * pow(Tk, 4) + (1 - C_ENV_EMISSION) * A;

      Q = C_BODY_ABSORPTION * G + A - E;
      switch (type) {
        case 1: { // sky
          B = Q > 0 ? -0.19 * Q : -0.32 * Q;
          Tk = Tk + (Q + B) / ((6.2 + 4.26 * wind) * (1 + 1 / C_BOWEN_RATIO));
          break;
        }
        case 3: { // building
          B = -0.55 * Q - 38;
          Tk = Tk + (Q + B) / ((6.2 + 4.26 * wind) * (1 + 1 / C_BOWEN_RATIO_BUILDING));
          break;
        }
        case 4: { //roads/impervious
          B = -0.70 * Q - 38;
          Tk = Tk + (Q + B) / ((6.2 + 4.26 * wind) * (1 + 1 / C_BOWEN_RATIO_IMPERVIOUS));
          break;
        }
        case 5: //pervious
        case 6: { //moving
          B = -0.34 * Q - 31;
          Tk = Tk + (Q + B) / ((6.2 + 4.26 * wind) * (1 + 1 / C_BOWEN_RATIO_PERVIOUS));
          break;
        }
      }
    }
    Ts = Tk - 273.15;
  }

	// VDI 3789 (25) TODO: Use
	double getReflectedGlobalRadiation(double globalRadiation) {
		return C_ALBEDO_ENV * globalRadiation;
	}

	double getCosMu(double alpha, double beta, double phi, double gamma) {
		double alphaR = alpha * M_PI / 180;
		double betaR = beta * M_PI / 180;
		double phiR = phi * M_PI / 180;
		double gammaR = gamma * M_PI / 180;
		return sin(gammaR)*cos(betaR) + cos(gammaR)*sin(betaR)*cos(alphaR - phiR);
	}

};
} // calculate

#endif // !CALCULATEMRT_MRTCALCULATOR_H_