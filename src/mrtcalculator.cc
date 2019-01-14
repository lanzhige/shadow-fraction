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

#ifdef _DEBUG
	debug_image_buffer = (unsigned char*)malloc(512 * 512 * 4 * sizeof(unsigned char));
	ilInit();
#endif // DEBUG
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
			//std::cout << shadow[((511 - i) * 512 + j) * 4]<<", "<< shadow[((511 - i) * 512 + j) * 4+1]<<", "<< shadow[((511 - i) * 512 + j) * 4+2] << std::endl;
			//float t = shadow[((511 - i) * 512 + j) * 4] + shadow[((511 - i) * 512 + j) * 4 + 1] + shadow[((511 - i) * 512 + j) * 4 + 2];
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
	@brief Initial the start accumulated temperature
	@param date The date info of input, including the day of year and the end hour of initializing
	@param data The fisheye data of the top direction. 
*/
void MRTCalculator::initTotalTemperature(Date &date, double lat, double lng, std::vector<unsigned char*> &seg, std::vector<unsigned char*> &depth) {
	for (int i = 0; i < 6; i++) {
		acc_tem_time[i].clear();
		acc_tem_time[i].resize(512 * 512, 0);
	}
	Date temp(date);
	while (temp < date) {
		double dZenithAngle, dAzimuth;
		// calculate sun position
		getSunPos(date, lat, lng, dZenithAngle, dAzimuth);
		std::vector<float*> *res = scSuf->calculateSurfaceImages(depth, dZenithAngle, dAzimuth);
		std::swap(res->at(0), res->at(3));
		std::swap(res->at(1), res->at(2));
		std::swap(res->at(2), res->at(4));
		std::swap(res->at(3), res->at(5));
		std::swap(res->at(4), res->at(5));

		for (int i = 0; i < 6; i++) {
			getSegShadow(seg[i], res->at(i), shaded_buffer[i], non_shaded_buffer[i], acc_tem_time[i]);
		}
		scSuf->deallocateSurfaceImages(res);
		temp += time_step;
	}
}

#ifdef _DEBUG
void MRTCalculator::saveImage(const unsigned char *data, std::string location, int size) {
	for (int i = 0; i < size; i++) {
		unsigned index = data[i];
		//std::cout << index << std::endl;
		for (int j = 0; j < 4; j++) {
			if (index > 6) debug_image_buffer[i * 4 + j] = 255;
			else debug_image_buffer[i * 4 + j] = color_map[index][j];
		}
	}

	ILuint imageID = ilGenImage();

	ilBindImage(imageID);
	ilTexImage(512, 512, 1, 4, IL_RGBA, IL_UNSIGNED_BYTE, debug_image_buffer);
	ilEnable(IL_FILE_OVERWRITE);
	ilSave(IL_JPG, &location[0]);
	ilDeleteImage(imageID);
}

void MRTCalculator::saveTemperatureImage(const vector<int> &data, std::string location) {
	for (int i = 0; i < data.size(); i++) {
		if (data[i]) {
			debug_image_buffer[i * 4] = 255 - data[i] / 4;
			debug_image_buffer[i * 4 + 1] = 0;
			debug_image_buffer[i * 4 + 2] = 0;
			debug_image_buffer[i * 4 + 3] = data[i] / 4;
		}
		else {
			for (int j = 0; j < 4; j++) debug_image_buffer[i * 4 + j] = 255;
		}
	}

	ILuint imageID = ilGenImage();

	ilBindImage(imageID);
	ilTexImage(512, 512, 1, 4, IL_RGBA, IL_UNSIGNED_BYTE, debug_image_buffer);
	ilEnable(IL_FILE_OVERWRITE);
	ilSave(IL_JPG, &location[0]);
	ilDeleteImage(imageID);
}
#endif //_DEBUG

double MRTCalculator::vectorAngle(double x, double y, double z, double azimuth, int direction) {
  double a_x, a_y, a_z;

  switch (direction) {
    case 0: {
      a_x = -cos(azimuth);
      a_y = 0;
      a_z = -sin(azimuth);
      break;
    }
    case 1: {
      a_x = cos(azimuth);
      a_y = 0;
      a_z = sin(azimuth);
      break;
    }
    case 2: {
      a_x = cos(azimuth);
      a_y = sin(azimuth);
      a_z = 0;
      break;
    }
    case 3: {
      a_x = cos(azimuth);
      a_y = -sin(azimuth);
      a_z = 0;
      break;
    }
    case 4: {
      a_x = sin(azimuth);
      a_y = 0;
      a_z = -cos(azimuth);
      break;
    }
    case 5: {
      a_x = -sin(azimuth);
      a_y = 0;
      a_z = cos(azimuth);
      break;
    }
  }

  double len_vect = sqrt(x * x + y * y + z * z);
  double mul = x * a_x + y * a_y + z * a_z;
  std::cout << "vector length: " << len_vect << " complement: " << mul << std::endl;
  return acos(mul / (len_vect));
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

	if (pre_day != date._day) {
		initTotalTemperature(date, geo.lat, geo.lng, seg_data, depth_data);
		pre_day = date._day;
	}

	// get sun exposed hours
	double sunHours = getSunHours(&date, data[0], geo.lat, geo.lng);

	double dZenithAngle = 0.0, dAzimuth = 0.0;
	// calculate sun position
	getSunPos(date, geo.lat, geo.lng, dZenithAngle, dAzimuth);
	// calculate if sun is visible
	bool sunVisible = getSunVisible(dZenithAngle, dAzimuth, data[0]);

	// shadow data
/*	std::cout << "zenith: " << dZenithAngle << " azimuth: " << dAzimuth << std::endl;
	
	double sum = 0, bias = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 512 * 512; j++) {
			sum += depth_data[i][j];
			if (j & 1) {
				bias -= depth_data[i][j];
			}
			else bias += depth_data[i][j];
		}
	}
	std::cout << "sum: " << sum << std::endl;
	std::cout << "bias: " << bias << std::endl;*/

	std::vector<float*> *res = scSuf->calculateSurfaceImages(depth_data, dZenithAngle, dAzimuth);

/*	sum = 0, bias = 0;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 512 * 512 * 4; j++) {
			sum += res->at(i)[j];
			if (j & 1) {
				bias -= res->at(i)[j];
			}
			else bias += res->at(i)[j];
		}
	}
	std::cout << "sum: " << sum << std::endl;
	std::cout << "bias: " << bias << std::endl;
	*/

	std::swap(res->at(0), res->at(3));
	std::swap(res->at(1), res->at(2));
	std::swap(res->at(2), res->at(4));
	std::swap(res->at(3), res->at(5));
	std::swap(res->at(4), res->at(5));

	for (int i = 0; i < 6; i++) {
		getSegShadow(seg_data[i], res->at(i), shaded_buffer[i], non_shaded_buffer[i], acc_tem_time[i]);
	}

	fisheye->getVF(shaded_buffer, shaded_fe_buffer, vfCalculator, shaded_vf);
	fisheye->getVF(non_shaded_buffer, non_shaded_fe_buffer, vfCalculator, non_shaded_vf);

#ifdef _DEBUG
	string name = date.to_string() + "_" + std::to_string(geo.lat) + "_" + std::to_string(geo.lng) + "_";
	for (int i = 0; i < 6; i++) {
		saveImage(shaded_fe_buffer[i], "./Debug_output/shaded_fisheye/" + name + debug_postfix[i], 512 * 512);
		saveImage(non_shaded_fe_buffer[i], "./Debug_output/nonshaded_fisheye/" + name + debug_postfix[i], 512 * 512);
		saveTemperatureImage(acc_tem_time[i], "./Debug_output/accumulated_temperature/" + name + debug_temp[i]);
	}
#endif // DEBUG

	double gamma = 90.0 - dZenithAngle;

	double phi = dAzimuth;

	double atmosphericPressure = 101325 * pow(1 - 2.25577 / 100000 * geo.altitude, 5.25588);

	// Relative Optical Air Mass
	double relativeOpticalAirMass = getRelativeOpticalAirMass(gamma);

	// Optical Thickness
	double opticalThickness = getOpticalThickness(relativeOpticalAirMass);
  if (gamma < 5) opticalThickness = 0.0548;
  if (gamma < 4) opticalThickness = 0.0519;
  if (gamma < 3) opticalThickness = 0.0491;
  if (gamma < 2) opticalThickness = 0.0463;
  if (gamma < 1) opticalThickness = 0.0435;
  if (gamma < 0) opticalThickness = 0.0408;

	// Irradiance of extra...
	double extraterrestialRadiation = getExtraterrestialRadiation(date.dayOfYear);
  std::cout << "I0: " << extraterrestialRadiation << std::endl;
  // Vapor pressure
  double vaporPressure = getVaporPressure(geo.relativeHumidity, geo.airTemperature);

  // Atmospheric Radiation
  double atmosphericRadiation = getAtmosphericRadiation(geo.airTemperature, vaporPressure, geo.cloudCover);

  // Direct Radiation
  double directRadiation = getDirectRadiation(extraterrestialRadiation, gamma, opticalThickness, relativeOpticalAirMass, atmosphericPressure, geo.cloudCover);

  std::cout << "IOrig: " << directRadiation <<std::endl;

	for (int direction = 0; direction < 6; direction++) {
		auto &surface = res->at(direction);
    for (int i = 0; i < 512; i++) {
      for (int j = 0; j < 512; j++) {
        std::cout << surface[i * 512 + j] << " " << surface[i * 512 + j + 1] << " " << surface[i * 512 + j + 2] << std::endl;
      }
    }
    std::cout << "direction: " << direction << std::endl;
		for (int i = 0; i < 512; i++) {
			for (int j = 0; j < 512; j++) {
				int shaded = static_cast<int>(surface[((511 - i) * 512 + j) * 4 + 3]);

				//Modelling radiation fluxes in simple and complex environments: basics of the RayMan model
				//--(4)
				double I;
				if (shaded == 1) {
					I = 0;
				}
				else {
          std::cout << "i: " << 511 - i << " j: " << j << " " << surface[((511 - i) * 512 + j) * 4 + 0] << " " << surface[((511 - i) * 512 + j) * 4 + 1] << " " << surface[((511 - i) * 512 + j) * 4 + 2] << " " << surface[((511 - i) * 512 + j) * 4 + 3] << std::endl;
          double xi = vectorAngle(surface[((511 - i) * 512 + j) * 4 + 0],
                                  surface[((511 - i) * 512 + j) * 4 + 1],
                                  surface[((511 - i) * 512 + j) * 4 + 2],
                                  dAzimuth,
                                  direction);
          xi = xi * 180 / pi;
          std::cout << "XI: " << xi <<std::endl;
          I = directRadiation * cos(gamma) / cos(dZenithAngle) * cos(xi);        
        }
        double G0 = getGlobalRadiation(extraterrestialRadiation, gamma, atmosphericPressure, 0);
        double directRadiationCloudless = getDirectRadiation(extraterrestialRadiation, gamma, opticalThickness, relativeOpticalAirMass, atmosphericPressure, 0);
        double transmittanceOfDirectRadiation = getTransmittanceOfDirectRadiation(directRadiationCloudless, extraterrestialRadiation);

        double D0 = getDiffuseRadiationCloudless(G0, directRadiationCloudless, transmittanceOfDirectRadiation, viewFactor[0], sunVisible);
        double D8 = getDiffuseRadiationOvercast(G0, viewFactor[0]);
        double D = getDiffuseRadiation(D0, D8, geo.cloudCover);
        std::cout << "I: " << I << " D: " << D <<std::endl;
        double G = I + D;
        double Ts, E;
        Ts = geo.airTemperature + 273.15;
        getSurfaceTemperature(Ts, E, atmosphericRadiation, G, geo.windVelocity, (int)seg_data[direction][i * 512 + j]);
        std::cout << "surface temperature: " << Ts <<std::endl;
        surf_temp[direction][i * 512 + j] = (unsigned char)Ts;
        surf_emission[direction][i * 512 + j] = (unsigned char)E;
        int pause;
        std::cin >> pause;
			}
		}
    
		/*std::cout << "direction: " << direction << std::endl;
		for (int i = 0; i < 512; i += 32) {
			for (int j = 0; j < 512; j += 32) {
				std::cout << "i: " << i << " j: " << j << " " << surface[i * 512 + j] << " " << surface[i * 512 + j + 1] << " " << surface[i * 512 + j + 2] << std::endl;
			}
		}*/
	}

  std::vector<std::vector<double>> useless = factors;
  fisheye->getVF(surf_temp, surf_temp_fisheye, vfCalculator, useless);
  fisheye->getVF(surf_emission, surf_emission_fisheye, vfCalculator, useless);

  std::vector<double> longwave(6, 0), diffuse(6, 0);
  vector<double> new_meanRT;
  for (int direction = 0; direction < 6; direction++) {
    std::vector<double> p(32, 0), p_e(32, 0);
    for (int i = 0; i < 512; i++) {
      for (int j = 0; j < 512; j++) {
        double d = sqrt((i + 0.5 - 256) * (i + 0.5 - 256) + (j + 0.5 - 256) * (j + 0.5 - 256));
        int ringIndex = (int) floor(d / 8.0);
        if (ringIndex < 32) {
          p[ringIndex] += surf_temp_fisheye[direction][i * 512 + j];
          p_e[ringIndex] += surf_emission_fisheye[direction][i * 512 + j];
        }
      }
    }
    double vf_t = 0, vf_e = 0;
    for (int i = 0; i < 32; i++) {
      vf_t += vfCalculator.svf_max[i] * p[i];
      vf_e += vfCalculator.svf_max[i] * p_e[i];
    }
    longwave[direction] = vf_t * pi / 64.0;
    diffuse[direction] = vf_e * pi / 64.0;
    new_meanRT.push_back(longwave[direction] + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuse[direction]);
  }

  double new_MRT_Sum = (0.22*new_meanRT[2] + 0.22*new_meanRT[3] + 0.22*new_meanRT[4] + 0.22 * new_meanRT[5] + 0.06 * new_meanRT[1] + 0.06 * new_meanRT[0]);

  double new_mrt = pow(new_MRT_Sum / C_STEFAN_BOLTZMANN, 0.25) - 273.15;
  std::cout << "new mrt: " << new_mrt << std::endl;
  exit(1);
	//------------------------------------------

	// Global Radiation
	double globalRadiation = getGlobalRadiation(extraterrestialRadiation, gamma, atmosphericPressure, geo.cloudCover);

	// Get Diffuse Radiation
	double globalRadiationCloudless = getGlobalRadiation(extraterrestialRadiation, gamma, atmosphericPressure, 0);
	double directRadiationCloudless = getDirectRadiation(extraterrestialRadiation, gamma, opticalThickness, relativeOpticalAirMass, atmosphericPressure, 0);
	double transmittanceOfDirectRadiation = getTransmittanceOfDirectRadiation(directRadiationCloudless, extraterrestialRadiation);

	double diffuseRadiationCloudless = getDiffuseRadiationCloudless(globalRadiationCloudless, directRadiationCloudless, transmittanceOfDirectRadiation, viewFactor[0], sunVisible);
	double diffuseRadiationOvercast = getDiffuseRadiationOvercast(globalRadiationCloudless, viewFactor[0]);
	double diffuseRadiation = getDiffuseRadiation(diffuseRadiationCloudless, diffuseRadiationOvercast, geo.cloudCover);



	// Iterative E and T_s
	std::pair<double, double> temp_Sky = getSurfaceTemperatureAndLongWaveRadiation(geo.airTemperature, atmosphericRadiation, globalRadiation, geo.windVelocity, "SKY");
	double surfaceTemperature_Sky = temp_Sky.first;
	double longWaveRadiation_Sky = temp_Sky.second;

	// Iterative E and T_s
	std::pair<double, double> temp_Tree = getSurfaceTemperatureAndLongWaveRadiation(geo.airTemperature, atmosphericRadiation, globalRadiation, geo.windVelocity, "TREE");
	double surfaceTemperature_Tree = temp_Tree.first;
	double longWaveRadiation_Tree = temp_Tree.second;

	// Iterative E and T_s
	std::pair<double, double> temp_Building = getSurfaceTemperatureAndLongWaveRadiation(geo.airTemperature, atmosphericRadiation, globalRadiation, geo.windVelocity, "BUILDING");
	double surfaceTemperature_Building = temp_Building.first;
	double longWaveRadiation_Building = temp_Building.second;

	// Iterative E and T_s
	std::pair<double, double> temp_Impervious = getSurfaceTemperatureAndLongWaveRadiation(geo.airTemperature, atmosphericRadiation, globalRadiation, geo.windVelocity, "IMPERVIOUS");
	double surfaceTemperature_Impervious = temp_Impervious.first;
	double longWaveRadiation_Impervious = temp_Impervious.second;

	// Iterative E and T_s
	std::pair<double, double> temp_Pervious = getSurfaceTemperatureAndLongWaveRadiation(geo.airTemperature, atmosphericRadiation, globalRadiation, geo.windVelocity, "PERVIOUS");
	double surfaceTemperature_Pervious = temp_Pervious.first;
	double longWaveRadiation_Pervious = temp_Pervious.second;

	// Iterative E and T_s
	std::pair<double, double> temp_MovingObjects = getSurfaceTemperatureAndLongWaveRadiation(geo.airTemperature, atmosphericRadiation, globalRadiation, geo.windVelocity, "MOVINGOBJECTS");
	double surfaceTemperature_MovingObjects = temp_MovingObjects.first;
	double longWaveRadiation_MovingObjects = temp_MovingObjects.second;

	std::vector<double> meanRT;
	std::vector<double> vf360(6, 0.0);
	std::vector<double> vfmrt(6, 0.0);
	for (auto& viewFactor : factors){
		double SurfaceTotal_Sky = viewFactor[0] * (atmosphericRadiation + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuseRadiation);
		double SurfaceTotal_Tree = viewFactor[1] * (longWaveRadiation_Tree + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuseRadiation);
		double SurfaceTotal_Building = viewFactor[2] * (longWaveRadiation_Building + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuseRadiation);
		double SurfaceTotal_Impervious = viewFactor[3] * (longWaveRadiation_Impervious + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuseRadiation);
		double SurfaceTotal_Pervious = viewFactor[4] * (longWaveRadiation_Pervious + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuseRadiation);
		double SurfaceTotal_MovingObjects = viewFactor[5] * (longWaveRadiation_MovingObjects + C_BODY_ABSORPTION / C_BODY_EMISSION * diffuseRadiation);

		for (int i = 0; i < 6; i++) {
			vf360[i] += viewFactor[i];
		}

		double AllRadiationFromFisheye = SurfaceTotal_Sky + SurfaceTotal_Tree + SurfaceTotal_Building + SurfaceTotal_Impervious + SurfaceTotal_Pervious + SurfaceTotal_MovingObjects;
		meanRT.push_back(AllRadiationFromFisheye);
	}

	double I_s = getDirectRadiationPerpendicular(extraterrestialRadiation, gamma, phi, opticalThickness, relativeOpticalAirMass, atmosphericPressure, geo.cloudCover);

	for (int i = 0; i < 6; i++) {
		vfmrt[i] = 0.22*(factors[2][i] + factors[3][i] + factors[4][i] + factors[5][i]) + (factors[0][i] + factors[1][i])*0.06;
	}

	double MRT_Sum = (0.22*meanRT[2] + 0.22*meanRT[3] + 0.22*meanRT[4] + 0.22*meanRT[5] + 0.06*meanRT[1] + 0.06*meanRT[0]);
	double g2 = gamma * gamma;
	double surfaceProjectionFactor = 0.000000381313131*g2*gamma - 0.000065530303030*g2 + 0.000302020202020*gamma + 0.308054545454531;

	if (sunVisible)
		MRT_Sum += C_BODY_ABSORPTION / C_BODY_EMISSION * surfaceProjectionFactor * I_s;

	double mrt = pow(MRT_Sum / C_STEFAN_BOLTZMANN, 0.25) - 273.15;

	scSuf->deallocateSurfaceImages(res);

	// generate output
	ss << geo.lat << "," << geo.lng << "," << date._year << "," << date._month << "," << date._day << "," << date._hour << "," << date.dayOfYear << "," << geo.airTemperature << "," << geo.relativeHumidity << "," << vaporPressure << "," << geo.windVelocity << ",";

	ss << sunHours << "," << sunVisible << "," << gamma << "," << phi << "," << atmosphericPressure << "," << relativeOpticalAirMass << "," <<
		opticalThickness << "," << extraterrestialRadiation << "," << directRadiation << "," << globalRadiation << "," <<
		diffuseRadiation << "," << atmosphericRadiation << ",";

	ss << surfaceTemperature_Sky << "," << surfaceTemperature_Tree << "," << surfaceTemperature_Building << "," << surfaceTemperature_Impervious << "," << surfaceTemperature_Pervious << "," << surfaceTemperature_MovingObjects << ",";
	ss << longWaveRadiation_Sky << "," << longWaveRadiation_Tree << "," << longWaveRadiation_Building << "," << longWaveRadiation_Impervious << "," << longWaveRadiation_Pervious << "," << longWaveRadiation_MovingObjects << ",";
	
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			ss << factors[i][j] << ",";
		}
	}
	for (auto&i : vf360) {
		ss << i/6.0 << ",";
	}

	for (auto&i : vfmrt) {
		ss << i << ",";
	}

	for (auto&i : meanRT) {
		ss << i << ",";
	}

	ss << mrt << ",";
	
	// calculate UTCI
	calculateUTCI(geo.airTemperature, geo.windVelocity, mrt, vaporPressure, ss);

	// calculate PET
	PETCalculator pet(geo.airTemperature, mrt, vaporPressure, geo.windVelocity, json);
	pet.calculate(&ss);

	for (int i = 0; i < 6; i++) {
		for (int j = 1; j < 6; j++) {
			ss << "," << shaded_vf[i][j];
		}
	}

	for (int i = 0; i < 6; i++) {
		for (int j = 1; j < 6; j++) {
			ss << ", " << non_shaded_vf[i][j];
		}
	}

	ss << "\n";
	*output += ss.str();
}
/**
	@brief Calculate the UTCI result
	@param Ta The air temperature in current time and position
	@param va The wind velocity in current time and positionggg
	@param Tmrt The mean radiant temperature
	@param ehPa The current vapor pressure
	@param ss The stringstream for output
*/
void MRTCalculator::calculateUTCI(double Ta, double va, double Tmrt, double ehPa, std::stringstream &ss) {
	double D_Tmrt = Tmrt - Ta;
	double Pa = ehPa / 10.0;//  convert vapour pressure to kPa

	double UTCI_approx = Ta +
		(0.607562052) +
		(-0.0227712343) * Ta +
		(8.06470249 * pow(10, (-4))) * Ta * Ta +
		(-1.54271372 * pow(10, (-4))) * Ta * Ta * Ta +
		(-3.24651735 * pow(10, (-6))) * Ta * Ta * Ta * Ta +
		(7.32602852 * pow(10, (-8))) * Ta * Ta * Ta * Ta * Ta +
		(1.35959073 * pow(10, (-9))) * Ta * Ta * Ta * Ta * Ta * Ta +
		(-2.25836520) * va +
		(0.0880326035) * Ta * va +
		(0.00216844454) * Ta * Ta * va +
		(-1.53347087 * pow(10, (-5))) * Ta * Ta * Ta * va +
		(-5.72983704 * pow(10, (-7))) * Ta * Ta * Ta * Ta * va +
		(-2.55090145 * pow(10, (-9))) * Ta * Ta * Ta * Ta * Ta * va +
		(-0.751269505) * va * va +
		(-0.00408350271) * Ta * va * va +
		(-5.21670675 * pow(10, (-5))) * Ta * Ta * va * va +
		(1.94544667 * pow(10, (-6))) * Ta * Ta * Ta * va * va +
		(1.14099531 * pow(10, (-8))) * Ta * Ta * Ta * Ta * va * va +
		(0.158137256) * va * va * va +
		(-6.57263143 * pow(10, (-5))) * Ta * va * va * va +
		(2.22697524 * pow(10, (-7))) * Ta * Ta * va * va * va +
		(-4.16117031 * pow(10, (-8))) * Ta * Ta * Ta * va * va * va +
		(-0.0127762753) * va * va * va * va +
		(9.66891875 * pow(10, (-6))) * Ta * va * va * va * va +
		(2.52785852 * pow(10, (-9))) * Ta * Ta * va * va * va * va +
		(4.56306672 * pow(10, (-4))) * va * va * va * va * va +
		(-1.74202546 * pow(10, (-7))) * Ta * va * va * va * va * va +
		(-5.91491269 * pow(10, (-6))) * va * va * va * va * va * va +
		(0.398374029) * D_Tmrt +
		(1.83945314 * pow(10, (-4))) * Ta * D_Tmrt +
		(-1.73754510 * pow(10, (-4))) * Ta * Ta * D_Tmrt +
		(-7.60781159 * pow(10, (-7))) * Ta * Ta * Ta * D_Tmrt +
		(3.77830287 * pow(10, (-8))) * Ta * Ta * Ta * Ta * D_Tmrt +
		(5.43079673 * pow(10, (-10))) * Ta * Ta * Ta * Ta * Ta * D_Tmrt +
		(-0.0200518269) * va * D_Tmrt +
		(8.92859837 * pow(10, (-4))) * Ta * va * D_Tmrt +
		(3.45433048 * pow(10, (-6))) * Ta * Ta * va * D_Tmrt +
		(-3.77925774 * pow(10, (-7))) * Ta * Ta * Ta * va * D_Tmrt +
		(-1.69699377 * pow(10, (-9))) * Ta * Ta * Ta * Ta * va * D_Tmrt +
		(1.69992415 * pow(10, (-4))) * va * va * D_Tmrt +
		(-4.99204314 * pow(10, (-5))) * Ta * va * va * D_Tmrt +
		(2.47417178 * pow(10, (-7))) * Ta * Ta * va * va * D_Tmrt +
		(1.07596466 * pow(10, (-8))) * Ta * Ta * Ta * va * va * D_Tmrt +
		(8.49242932 * pow(10, (-5))) * va * va * va * D_Tmrt +
		(1.35191328 * pow(10, (-6))) * Ta * va * va * va * D_Tmrt +
		(-6.21531254 * pow(10, (-9))) * Ta * Ta * va * va * va * D_Tmrt +
		(-4.99410301 * pow(10, (-6))) * va * va * va * va * D_Tmrt +
		(-1.89489258 * pow(10, (-8))) * Ta * va * va * va * va * D_Tmrt +
		(8.15300114 * pow(10, (-8))) * va * va * va * va * va * D_Tmrt +
		(7.55043090 * pow(10, (-4))) * D_Tmrt * D_Tmrt +
		(-5.65095215 * pow(10, (-5))) * Ta * D_Tmrt * D_Tmrt +
		(-4.52166564 * pow(10, (-7))) * Ta * Ta * D_Tmrt * D_Tmrt +
		(2.46688878 * pow(10, (-8))) * Ta * Ta * Ta * D_Tmrt * D_Tmrt +
		(2.42674348 * pow(10, (-10))) * Ta * Ta * Ta * Ta * D_Tmrt * D_Tmrt +
		(1.54547250 * pow(10, (-4))) * va * D_Tmrt * D_Tmrt +
		(5.24110970 * pow(10, (-6))) * Ta * va * D_Tmrt * D_Tmrt +
		(-8.75874982 * pow(10, (-8))) * Ta * Ta * va * D_Tmrt * D_Tmrt +
		(-1.50743064 * pow(10, (-9))) * Ta * Ta * Ta * va * D_Tmrt * D_Tmrt +
		(-1.56236307 * pow(10, (-5))) * va * va * D_Tmrt * D_Tmrt +
		(-1.33895614 * pow(10, (-7))) * Ta * va * va * D_Tmrt * D_Tmrt +
		(2.49709824 * pow(10, (-9))) * Ta * Ta * va * va * D_Tmrt * D_Tmrt +
		(6.51711721 * pow(10, (-7))) * va * va * va * D_Tmrt * D_Tmrt +
		(1.94960053 * pow(10, (-9))) * Ta * va * va * va * D_Tmrt * D_Tmrt +
		(-1.00361113 * pow(10, (-8))) * va * va * va * va * D_Tmrt * D_Tmrt +
		(-1.21206673 * pow(10, (-5))) * D_Tmrt * D_Tmrt * D_Tmrt +
		(-2.18203660 * pow(10, (-7))) * Ta * D_Tmrt * D_Tmrt * D_Tmrt +
		(7.51269482 * pow(10, (-9))) * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt +
		(9.79063848 * pow(10, (-11))) * Ta * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt +
		(1.25006734 * pow(10, (-6))) * va * D_Tmrt * D_Tmrt * D_Tmrt +
		(-1.81584736 * pow(10, (-9))) * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt +
		(-3.52197671 * pow(10, (-10))) * Ta * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt +
		(-3.36514630 * pow(10, (-8))) * va * va * D_Tmrt * D_Tmrt * D_Tmrt +
		(1.35908359 * pow(10, (-10))) * Ta * va * va * D_Tmrt * D_Tmrt * D_Tmrt +
		(4.17032620 * pow(10, (-10))) * va * va * va * D_Tmrt * D_Tmrt * D_Tmrt +
		(-1.30369025 * pow(10, (-9))) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(4.13908461 * pow(10, (-10))) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(9.22652254 * pow(10, (-12))) * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(-5.08220384 * pow(10, (-9))) * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(-2.24730961 * pow(10, (-11))) * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(1.17139133 * pow(10, (-10))) * va * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(6.62154879 * pow(10, (-10))) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(4.03863260 * pow(10, (-13))) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(1.95087203 * pow(10, (-12))) * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(-4.73602469 * pow(10, (-12))) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt +
		(5.12733497) * Pa +
		(-0.312788561) * Ta * Pa +
		(-0.0196701861) * Ta * Ta * Pa +
		(9.99690870 * pow(10, (-4))) * Ta * Ta * Ta * Pa +
		(9.51738512 * pow(10, (-6))) * Ta * Ta * Ta * Ta * Pa +
		(-4.66426341 * pow(10, (-7))) * Ta * Ta * Ta * Ta * Ta * Pa +
		(0.548050612) * va * Pa +
		(-0.00330552823) * Ta * va * Pa +
		(-0.00164119440) * Ta * Ta * va * Pa +
		(-5.16670694 * pow(10, (-6))) * Ta * Ta * Ta * va * Pa +
		(9.52692432 * pow(10, (-7))) * Ta * Ta * Ta * Ta * va * Pa +
		(-0.0429223622) * va * va * Pa +
		(0.00500845667) * Ta * va * va * Pa +
		(1.00601257 * pow(10, (-6))) * Ta * Ta * va * va * Pa +
		(-1.81748644 * pow(10, (-6))) * Ta * Ta * Ta * va * va * Pa +
		(-1.25813502 * pow(10, (-3))) * va * va * va * Pa +
		(-1.79330391 * pow(10, (-4))) * Ta * va * va * va * Pa +
		(2.34994441 * pow(10, (-6))) * Ta * Ta * va * va * va * Pa +
		(1.29735808 * pow(10, (-4))) * va * va * va * va * Pa +
		(1.29064870 * pow(10, (-6))) * Ta * va * va * va * va * Pa +
		(-2.28558686 * pow(10, (-6))) * va * va * va * va * va * Pa +
		(-0.0369476348) * D_Tmrt * Pa +
		(0.00162325322) * Ta * D_Tmrt * Pa +
		(-3.14279680 * pow(10, (-5))) * Ta * Ta * D_Tmrt * Pa +
		(2.59835559 * pow(10, (-6))) * Ta * Ta * Ta * D_Tmrt * Pa +
		(-4.77136523 * pow(10, (-8))) * Ta * Ta * Ta * Ta * D_Tmrt * Pa +
		(8.64203390 * pow(10, (-3))) * va * D_Tmrt * Pa +
		(-6.87405181 * pow(10, (-4))) * Ta * va * D_Tmrt * Pa +
		(-9.13863872 * pow(10, (-6))) * Ta * Ta * va * D_Tmrt * Pa +
		(5.15916806 * pow(10, (-7))) * Ta * Ta * Ta * va * D_Tmrt * Pa +
		(-3.59217476 * pow(10, (-5))) * va * va * D_Tmrt * Pa +
		(3.28696511 * pow(10, (-5))) * Ta * va * va * D_Tmrt * Pa +
		(-7.10542454 * pow(10, (-7))) * Ta * Ta * va * va * D_Tmrt * Pa +
		(-1.24382300 * pow(10, (-5))) * va * va * va * D_Tmrt * Pa +
		(-7.38584400 * pow(10, (-9))) * Ta * va * va * va * D_Tmrt * Pa +
		(2.20609296 * pow(10, (-7))) * va * va * va * va * D_Tmrt * Pa +
		(-7.32469180 * pow(10, (-4))) * D_Tmrt * D_Tmrt * Pa +
		(-1.87381964 * pow(10, (-5))) * Ta * D_Tmrt * D_Tmrt * Pa +
		(4.80925239 * pow(10, (-6))) * Ta * Ta * D_Tmrt * D_Tmrt * Pa +
		(-8.75492040 * pow(10, (-8))) * Ta * Ta * Ta * D_Tmrt * D_Tmrt * Pa +
		(2.77862930 * pow(10, (-5))) * va * D_Tmrt * D_Tmrt * Pa +
		(-5.06004592 * pow(10, (-6))) * Ta * va * D_Tmrt * D_Tmrt * Pa +
		(1.14325367 * pow(10, (-7))) * Ta * Ta * va * D_Tmrt * D_Tmrt * Pa +
		(2.53016723 * pow(10, (-6))) * va * va * D_Tmrt * D_Tmrt * Pa +
		(-1.72857035 * pow(10, (-8))) * Ta * va * va * D_Tmrt * D_Tmrt * Pa +
		(-3.95079398 * pow(10, (-8))) * va * va * va * D_Tmrt * D_Tmrt * Pa +
		(-3.59413173 * pow(10, (-7))) * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(7.04388046 * pow(10, (-7))) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(-1.89309167 * pow(10, (-8))) * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(-4.79768731 * pow(10, (-7))) * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(7.96079978 * pow(10, (-9))) * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(1.62897058 * pow(10, (-9))) * va * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(3.94367674 * pow(10, (-8))) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(-1.18566247 * pow(10, (-9))) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(3.34678041 * pow(10, (-10))) * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(-1.15606447 * pow(10, (-10))) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa +
		(-2.80626406) * Pa * Pa +
		(0.548712484) * Ta * Pa * Pa +
		(-0.00399428410) * Ta * Ta * Pa * Pa +
		(-9.54009191 * pow(10, (-4))) * Ta * Ta * Ta * Pa * Pa +
		(1.93090978 * pow(10, (-5))) * Ta * Ta * Ta * Ta * Pa * Pa +
		(-0.308806365) * va * Pa * Pa +
		(0.0116952364) * Ta * va * Pa * Pa +
		(4.95271903 * pow(10, (-4))) * Ta * Ta * va * Pa * Pa +
		(-1.90710882 * pow(10, (-5))) * Ta * Ta * Ta * va * Pa * Pa +
		(0.00210787756) * va * va * Pa * Pa +
		(-6.98445738 * pow(10, (-4))) * Ta * va * va * Pa * Pa +
		(2.30109073 * pow(10, (-5))) * Ta * Ta * va * va * Pa * Pa +
		(4.17856590 * pow(10, (-4))) * va * va * va * Pa * Pa +
		(-1.27043871 * pow(10, (-5))) * Ta * va * va * va * Pa * Pa +
		(-3.04620472 * pow(10, (-6))) * va * va * va * va * Pa * Pa +
		(0.0514507424) * D_Tmrt * Pa * Pa +
		(-0.00432510997) * Ta * D_Tmrt * Pa * Pa +
		(8.99281156 * pow(10, (-5))) * Ta * Ta * D_Tmrt * Pa * Pa +
		(-7.14663943 * pow(10, (-7))) * Ta * Ta * Ta * D_Tmrt * Pa * Pa +
		(-2.66016305 * pow(10, (-4))) * va * D_Tmrt * Pa * Pa +
		(2.63789586 * pow(10, (-4))) * Ta * va * D_Tmrt * Pa * Pa +
		(-7.01199003 * pow(10, (-6))) * Ta * Ta * va * D_Tmrt * Pa * Pa +
		(-1.06823306 * pow(10, (-4))) * va * va * D_Tmrt * Pa * Pa +
		(3.61341136 * pow(10, (-6))) * Ta * va * va * D_Tmrt * Pa * Pa +
		(2.29748967 * pow(10, (-7))) * va * va * va * D_Tmrt * Pa * Pa +
		(3.04788893 * pow(10, (-4))) * D_Tmrt * D_Tmrt * Pa * Pa +
		(-6.42070836 * pow(10, (-5))) * Ta * D_Tmrt * D_Tmrt * Pa * Pa +
		(1.16257971 * pow(10, (-6))) * Ta * Ta * D_Tmrt * D_Tmrt * Pa * Pa +
		(7.68023384 * pow(10, (-6))) * va * D_Tmrt * D_Tmrt * Pa * Pa +
		(-5.47446896 * pow(10, (-7))) * Ta * va * D_Tmrt * D_Tmrt * Pa * Pa +
		(-3.59937910 * pow(10, (-8))) * va * va * D_Tmrt * D_Tmrt * Pa * Pa +
		(-4.36497725 * pow(10, (-6))) * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa +
		(1.68737969 * pow(10, (-7))) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa +
		(2.67489271 * pow(10, (-8))) * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa +
		(3.23926897 * pow(10, (-9))) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa +
		(-0.0353874123) * Pa * Pa * Pa +
		(-0.221201190) * Ta * Pa * Pa * Pa +
		(0.0155126038) * Ta * Ta * Pa * Pa * Pa +
		(-2.63917279 * pow(10, (-4))) * Ta * Ta * Ta * Pa * Pa * Pa +
		(0.0453433455) * va * Pa * Pa * Pa +
		(-0.00432943862) * Ta * va * Pa * Pa * Pa +
		(1.45389826 * pow(10, (-4))) * Ta * Ta * va * Pa * Pa * Pa +
		(2.17508610 * pow(10, (-4))) * va * va * Pa * Pa * Pa +
		(-6.66724702 * pow(10, (-5))) * Ta * va * va * Pa * Pa * Pa +
		(3.33217140 * pow(10, (-5))) * va * va * va * Pa * Pa * Pa +
		(-0.00226921615) * D_Tmrt * Pa * Pa * Pa +
		(3.80261982 * pow(10, (-4))) * Ta * D_Tmrt * Pa * Pa * Pa +
		(-5.45314314 * pow(10, (-9))) * Ta * Ta * D_Tmrt * Pa * Pa * Pa +
		(-7.96355448 * pow(10, (-4))) * va * D_Tmrt * Pa * Pa * Pa +
		(2.53458034 * pow(10, (-5))) * Ta * va * D_Tmrt * Pa * Pa * Pa +
		(-6.31223658 * pow(10, (-6))) * va * va * D_Tmrt * Pa * Pa * Pa +
		(3.02122035 * pow(10, (-4))) * D_Tmrt * D_Tmrt * Pa * Pa * Pa +
		(-4.77403547 * pow(10, (-6))) * Ta * D_Tmrt * D_Tmrt * Pa * Pa * Pa +
		(1.73825715 * pow(10, (-6))) * va * D_Tmrt * D_Tmrt * Pa * Pa * Pa +
		(-4.09087898 * pow(10, (-7))) * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa * Pa +
		(0.614155345) * Pa * Pa * Pa * Pa +
		(-0.0616755931) * Ta * Pa * Pa * Pa * Pa +
		(0.00133374846) * Ta * Ta * Pa * Pa * Pa * Pa +
		(0.00355375387) * va * Pa * Pa * Pa * Pa +
		(-5.13027851 * pow(10, (-4))) * Ta * va * Pa * Pa * Pa * Pa +
		(1.02449757 * pow(10, (-4))) * va * va * Pa * Pa * Pa * Pa +
		(-0.00148526421) * D_Tmrt * Pa * Pa * Pa * Pa +
		(-4.11469183 * pow(10, (-5))) * Ta * D_Tmrt * Pa * Pa * Pa * Pa +
		(-6.80434415 * pow(10, (-6))) * va * D_Tmrt * Pa * Pa * Pa * Pa +
		(-9.77675906 * pow(10, (-6))) * D_Tmrt * D_Tmrt * Pa * Pa * Pa * Pa +
		(0.0882773108) * Pa * Pa * Pa * Pa * Pa +
		(-0.00301859306) * Ta * Pa * Pa * Pa * Pa * Pa +
		(0.00104452989) * va * Pa * Pa * Pa * Pa * Pa +
		(2.47090539 * pow(10, (-4))) * D_Tmrt * Pa * Pa * Pa * Pa * Pa +
		(0.00148348065) * Pa * Pa * Pa * Pa * Pa * Pa;
		ss << UTCI_approx << ",";
}

} // calculate