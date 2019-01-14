#include "jsonstruct.h"

#include<vector>
#include<iostream>

namespace file {
JsonStruct::JsonStruct(){
	altitude = new std::vector<double>();
	airTemperature = new std::vector<double>();
	relativeHumidity = new std::vector<double>();
	windVelocity = new std::vector<double>();
	cloudCover = new std::vector<double>();

	tiles = new std::vector<std::string>();
}

JsonStruct::~JsonStruct() {
	if (altitude) delete altitude;
	if (airTemperature) delete airTemperature;
	if (relativeHumidity) delete relativeHumidity;
	if (windVelocity) delete windVelocity;
	if (cloudCover) delete cloudCover;

	if (tiles) delete tiles;
}
}