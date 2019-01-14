#include <iostream>
#include <string>

#include "fileio.h"
#include "jsonstruct.h"
#include "jsonreader.h"

int main(int argc, char** argv) {
	std::clock_t start;
	start = std::clock();
	std::cout << "start openMRT program" << std::endl;
	file::JsonStruct json;
	std::string path("./config.json");
	file::JsonReader reader;
	reader.JsonParser(path, json);
	file::FileIO f(&json);
	f.load();
	double duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
	std::cout << "seconds used: " << duration << std::endl;
	return 0;
}