#include "mrtcalculator.h"

#include <sstream>

#include"IL/il.h"

#include "petcalculator.h"

namespace calculate {
/**
	@brief Initialize shadow calculator
*/
void MRTCalculator::init(int timeStep) {
	if (scSuf == nullptr) {
		std::cout << "initializing surface calculator" << std::endl;
		scSuf = new SC_SurfaceCalc();
		scSuf->parameters.WINDOW_INITIAL_HEIGHT = 512;
		scSuf->parameters.WINDOW_INITIAL_WIDTH = 512;
		scSuf->parameters.SHADOW_MAP_EXPORT = false;
		scSuf->parameters.DISPLAY_OUTPUT_IN_WINDOW = false;
		scSuf->initialize();
	}

	for (auto p : shaded_buffer) {
		if (p != nullptr) free(p);
	}
	//shaded_buffer.clear();
	for (auto p : non_shaded_buffer){
		if (p != nullptr) free(p);
	}

	for (auto p : shaded_fe_buffer) {
		if (p != nullptr) free(p);
	}

	//shaded_buffer.clear();
	for (auto p : non_shaded_fe_buffer) {
		if (p != nullptr) free(p);
	}

	//non_shaded_buffer.clear();
	shaded_buffer.resize(6);
	non_shaded_buffer.resize(6);
	shaded_fe_buffer.resize(6);
	non_shaded_fe_buffer.resize(6);
	for (int i = 0; i < 6; i++) {
		shaded_buffer[i] = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));
		non_shaded_buffer[i] = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));
		shaded_fe_buffer[i] = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));
		non_shaded_fe_buffer[i] = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));
	}
	
	acc_tem_time.resize(6);
	for (int i = 0; i < 6; i++) {
		shaded_vf.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
		non_shaded_vf.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
		acc_tem_time[i].resize(512 * 512, 0);
	}

	time_step = timeStep;
	int pre_day = -1;

  surf_temp.resize(6);
  surf_emission.resize(6);
  surf_temp_fisheye.resize(6);
  surf_emission_fisheye.resize(6);
  for (int i = 0; i < 6; i++) {
    surf_temp[i] = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));
    surf_emission[i] = (unsigned char*)malloc(512 * 512 * sizeof(unsigned char));
    surf_temp_fisheye[i] = (unsigned char*)malloc(512 * 512 * sizeof(unsigned char));
    surf_emission_fisheye[i] = (unsigned char*)malloc(512 * 512 * sizeof(unsigned char));
  }
  fc = new FractionCalculator(512, 512);
}

/**
	@brief Get the hours exposed to sun
	@param date The date object which contains date info of current time stamp
	@param data The segmented binary data buffer
	@param lat The latitude
	@param lng The longitude
*/
double MRTCalculator::getSunHours(Date *date, const unsigned char *data, double lat, double lng) {
	// use temporary variable to store the original date info
	int tempHour = date->_hour, tempMinute = date->_minute;
	date->sunHours = 0;
	for (date->_hour = 0; date->_hour<24; (date->_hour)++) {
		for (date->_minute = 0; date->_minute<60; date->_minute += 12) {
			double dZenithAngle, dAzimuth;
			getSunPos(*date, lat, lng, dZenithAngle, dAzimuth);
			if (getSunVisible(dZenithAngle, dAzimuth, data))
				date->sunHours += 0.2;
		}
	}
	date->_hour = tempHour;
	date->_minute = tempMinute;

	return date->sunHours;
}

/**
	@brief Get the current sun position
	@param date The date object which contains date info of current time stamp
	@param lat The latitude
	@param lng The longitude
	@param dZenithAngle The Zenith Angle of sun position as output
	@param dAzimuth The Azimuth angle of sun position as output
*/
void MRTCalculator::getSunPos(const Date &date, double lat, double lng, double &dZenithAngle, double &dAzimuth) {
	// Calculate time of the day in UT decimal hours
	double dDecimalHours = date._hour + date._minute / 60.0;
	// Calculate current Julian Day
	double liAux1 = (date._month - 14) / 12;
	double liAux2 = (1461.0 * (date._year + 4800.0 + liAux1)) / 4.0 + (367.0 * (date._month
		- 2.0 - 12.0 * liAux1)) / 12.0 - (3.0 * ((date._year + 4900
			+ liAux1) / 100.0)) / 4.0 + date._day - 32075.0;
	double dJulianDate = liAux2 - 0.5 + dDecimalHours / 24.0;
	// Calculate difference between current Julian Day and JD 2451545.0
	double dElapsedJulianDays = dJulianDate - 2451545.0;

	// Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
	// ecliptic in radians but without limiting the angle to be less than 2*Pi
	// (i.e., the result may be greater than 2*Pi)
	double dOmega = 2.1429 - 0.0010394594*dElapsedJulianDays;
	double dMeanLongitude = 4.8950630 + 0.017202791698*dElapsedJulianDays; // Radians
	double dMeanAnomaly = 6.2400600 + 0.0172019699*dElapsedJulianDays;
	double dEclipticLongitude = dMeanLongitude + 0.03341607*sin(dMeanAnomaly)
		+ 0.00034894*sin(2.0 * dMeanAnomaly) - 0.0001134
		- 0.0000203*sin(dOmega);
	double dEclipticObliquity = 0.4090928 - 6.2140e-9*dElapsedJulianDays
		+ 0.0000396*cos(dOmega);

	// Calculate celestial coordinates ( right ascension and declination ) in radians
	// but without limiting the angle to be less than 2*Pi (i.e., the result may be
	// greater than 2*Pi)
	double dSin_EclipticLongitude = sin(dEclipticLongitude);
	double dY = cos(dEclipticObliquity) * dSin_EclipticLongitude;
	double dX = cos(dEclipticLongitude);
	double dRightAscension = atan2(dY, dX);
	if (dRightAscension < 0.0) dRightAscension = dRightAscension + twopi;
	double dDeclination = asin(sin(dEclipticObliquity)*dSin_EclipticLongitude);

	// Calculate local coordinates ( azimuth and zenith angle ) in degrees
	double dGreenwichMeanSiderealTime = 6.6974243242 +
		0.0657098283*dElapsedJulianDays
		+ dDecimalHours;
	double dLocalMeanSiderealTime = (dGreenwichMeanSiderealTime * 15
		+ lng) * rad;
	double dHourAngle = dLocalMeanSiderealTime - dRightAscension;
	double dLatitudeInRadians = lat*rad;
	double dCos_Latitude = cos(dLatitudeInRadians);
	double dSin_Latitude = sin(dLatitudeInRadians);
	double dCos_HourAngle = cos(dHourAngle);

	dZenithAngle = (acos(dCos_Latitude*dCos_HourAngle
		*cos(dDeclination) + sin(dDeclination)*dSin_Latitude));

	dY = -sin(dHourAngle);
	dX = tan(dDeclination) * dCos_Latitude - dSin_Latitude * dCos_HourAngle;
	dAzimuth = atan2(dY, dX);
	if (dAzimuth < 0.0)	dAzimuth = dAzimuth + twopi;
	dAzimuth = dAzimuth / rad;
	// Parallax Correction
	double dParallax = (dEarthMeanRadius / dAstronomicalUnit) * sin(dZenithAngle);
	dZenithAngle = (dZenithAngle + dParallax) / rad;
}

/**
	@brief Transform Zenith angle and Azimuth angle to pixel coordinates
	@param dZenithAngle The Zenith Angle of sun position
	@param dAzimuth The Azimuth angle of sun position
	@param i The pixel coordinate of x axis as output
	@param j The pixel coordinate of y axis as output
*/
void MRTCalculator::sunPositionToImageCoordinates(double dZenithAngle, double dAzimuth, int &i, int &j) {
	double zaRad = dZenithAngle / 180 * M_PI;
	double aRad = dAzimuth / 180 * M_PI;
	double x = sin(zaRad)*cos(aRad);
	double y = sin(zaRad)*sin(aRad);
	double z = cos(zaRad);

	double r = atan2(sqrt(x*x + y * y), z) / M_PI;
	double phi = atan2(y, x);

	i = static_cast<int>(512 - (r * cos(phi) + 0.5) * 512);
	j = static_cast<int>((r * sin(phi) + 0.5) * 512);
}

/**
	@brief Judge if the sun is visible
	@param dZenithAngle The Zenith Angle of sun position
	@param dAzimuth The Azimuth angle of sun position
	@param data Segmented binary data buffer
	@return Returns boolean indicates sun visible or not
*/
bool MRTCalculator::getSunVisible(double dZenithAngle, double dAzimuth, const unsigned char *data) {
	int black = 0;
	int white = 0;
	if (dZenithAngle >= 0 && dZenithAngle <= 90) {
		int i, j;
		sunPositionToImageCoordinates(dZenithAngle, dAzimuth, i, j);
		int index_i;
		if (i<256) {
			index_i = 256 + (256 - i);
		}
		else if (i>256) {
			index_i = 256 - (i - 256);
		}
		else
			index_i = 256;

		unsigned char pixels[] = { 
			data[512 * (index_i - 1) + (j - 1)], 
			data[512 * (index_i - 1) + (j)],
			data[512 * (index_i - 1) + (j + 1)],
			data[512 * (index_i)+(j - 1)], 
			data[512 * (index_i)+(j)],
			data[512 * (index_i)+(j + 1)],
			data[512 * (index_i + 1) + (j + 1)],
			data[512 * (index_i + 1) + (j)],
			data[512 * (index_i + 1) + (j + 1)] };

		for (int i = 0; i<9; i++) {
			if (pixels[i] == 1) {
				white++;
			}
			else {
				black++;
			}
		}
	}
	return (white > black);
}

void MRTCalculator::getSegShadow(unsigned char* segmented, float* shadow, unsigned char* shaded, unsigned char* non_shaded, vector<int> &time) {
	for (int i = 0; i < 512; i++) {
		for (int j = 0; j < 512; j++) {
			int temp = static_cast<int>(shadow[((511 - i) * 512 + j) * 4 + 3]);
			if (temp && (int)segmented[i * 512 + j] != 1) {
				shaded[i * 512 + j] = segmented[i * 512 + j];
				non_shaded[i * 512 + j] = 0;
			}
			else {
				non_shaded[i * 512 + j] = segmented[i * 512 + j];
				time[(511 - i) * 512 + j] += time_step;
				shaded[i * 512 + j] = 0;
			}
		}
	}
}

/**
	@brief Calculate the Mean Radiant Temperature
	@param date The date object which contains date info of current time stamp
	@param data The fisheye data buffer
	@param geo The geo structure contains latitude and longitude information
	@param viewFactor The sky view factor
	@param output The output buffer
	@param factors The view factors for six direction
	@param json The json structure containing setting configurations
	@param depth_data The depth data buffer
*/
void MRTCalculator::calculate(Date &date, const std::vector<unsigned char *> &data, const file::GeoStruct &geo, const std::vector<double> &viewFactor, std::string *output, const std::vector<std::vector<double>> &factors, const file::JsonStruct *json, std::vector<unsigned char*>& depth_data, Fisheye *fisheye, VFCalculator &vfCalculator, std::vector<unsigned char*>& seg_data) {
	std::stringstream ss;
	ss.precision(13);
	// get sun exposed hours
	double sunHours = getSunHours(&date, data[0], geo.lat, geo.lng);

	double dZenithAngle = 0.0, dAzimuth = 0.0;
	// calculate sun position
	getSunPos(date, geo.lat, geo.lng, dZenithAngle, dAzimuth);
	// calculate if sun is visible
	bool sunVisible = getSunVisible(dZenithAngle, dAzimuth, data[0]);

	std::vector<float*> *res = scSuf->calculateSurfaceImages(depth_data, dZenithAngle, dAzimuth);

	std::swap(res->at(0), res->at(3));
	std::swap(res->at(1), res->at(2));
	std::swap(res->at(2), res->at(4));
	std::swap(res->at(3), res->at(5));
	std::swap(res->at(4), res->at(5));

  ss << geo.lat << "," << geo.lng << "," << date._year << "," << date._month << "," << date._day << "," << date._hour << "," << date.dayOfYear << ",";
  ss << sunHours << "," << (int)sunVisible << ",";

  double sum = 0;
	for (int i = 0; i < 6; i++) {
    std::vector<double> fraction = fc->getFraction(seg_data[i], res->at(i));
    std::cout<< "direction: "<< i <<std::endl;
    for (int j = 1; j < fraction.size(); j++) {
      std::cout << fraction[j] << " ";
      
      if (j == 1) {
        sum += fraction[j];
        sum += fraction[7];
        ss << fraction[j] << ",";
      } else if (j != 7) {
        sum += fraction[j];
        ss << fraction[j] << ",";
      }
    }
    std::cout<< " sum :" << sum << std::endl;
	}

	scSuf->deallocateSurfaceImages(res);
	
	ss << "\n";
	*output += ss.str();
}

} // calculate