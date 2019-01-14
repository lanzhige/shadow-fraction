#ifndef CALCULATEMRT_JSONSTRUCT_H_
#define CALCULATEMRT_JSONSTRUCT_H_

#include <string>
#include <vector>

namespace file {
struct JsonStruct {
	/*path*/
	std::string iRoot;
	std::string oRoot;

	/*time*/
	std::string startTime;
	std::string endTime;
	int timeStep;

	/*measurements*/
	std::vector<double> *altitude;
	std::vector<double> *airTemperature;
	std::vector<double> *relativeHumidity;
	std::vector<double> *windVelocity;
	std::vector<double> *cloudCover;

	/*subjects*/
	double height;
	double weight;
	double age;
	int male;
	double clothing;
	double activity;
	int standing;
	
	/*tile folders*/
	std::vector<std::string> *tiles;
	JsonStruct();
	~JsonStruct();
};
}

#endif // !CALCULATEMRT_JSONSTRUCT_H_