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