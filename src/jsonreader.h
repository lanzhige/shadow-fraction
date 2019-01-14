#ifndef CALCULATEMRT_JSONREADER_H_
#define CALCULATEMRT_JSONREADER_H_

#include <string>
#include <vector>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include "jsonstruct.h"

namespace file {
class JsonReader {
 public:
  JsonReader();
  void JsonParser(const std::string &path, JsonStruct &jsonStruct);
  ~JsonReader();

 private:
  template <typename T>
  void configEntryToVector(boost::property_tree::ptree
	  const& pt,
	  boost::property_tree::ptree::key_type const& key,
	  std::vector<T> *data);
};
} // namespace file

#endif // !CALCULATEMRT_JSONREADER_H_
