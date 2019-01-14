#include "jsonreader.h"

#include <string>
#include <vector>
#include <iostream>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

namespace file {
JsonReader::JsonReader() {
}

void JsonReader::JsonParser(const std::string &path, JsonStruct &jsonStruct) {
  std::ifstream jsonFile(path);
  std::cout << "loading json configurations" << std::endl;
  boost::property_tree::ptree pt;
  boost::property_tree::read_json(jsonFile, pt);

  jsonStruct.iRoot = pt.get<std::string>("iRoot");
  jsonStruct.oRoot = pt.get<std::string>("oRoot");

  auto time = pt.get_child("time").begin();
  jsonStruct.startTime = time->second.get_value<std::string>();
  time++;
  jsonStruct.endTime = time->second.get_value<std::string>();
  time++;
  jsonStruct.timeStep = time->second.get_value<int>();
  /*
  auto measurements = pt.get_child("measurements");
  configEntryToVector<double>(pt, "measurements.altitude", jsonStruct.altitude);
  configEntryToVector<double>(pt, "measurements.airTemperature", jsonStruct.airTemperature);
  configEntryToVector<double>(pt, "measurements.relativeHumidity", jsonStruct.relativeHumidity);
  configEntryToVector<double>(pt, "measurements.windVelocity", jsonStruct.windVelocity);
  configEntryToVector<double>(pt, "measurements.cloudCover", jsonStruct.cloudCover);

  auto subject = pt.get_child("subject");
  jsonStruct.height = subject.get<double>("height");
  jsonStruct.weight = subject.get<double>("weight");
  jsonStruct.age = subject.get<double>("age");
  jsonStruct.male = subject.get<int>("male");
  jsonStruct.clothing = subject.get<double>("clothing");
  jsonStruct.activity = subject.get<double>("activity");
  jsonStruct.standing = subject.get<int>("standing");

  configEntryToVector<std::string>(pt, "tiles", jsonStruct.tiles);
  */
}

template <typename T>
void JsonReader::configEntryToVector(boost::property_tree::ptree
    const& pt,
    boost::property_tree::ptree::key_type const& key,
	std::vector<T> *data) {
  for (auto& item : pt.get_child(key))
    data->push_back(item.second.get_value<T>());
}

JsonReader::~JsonReader() {
}
} // namespace file