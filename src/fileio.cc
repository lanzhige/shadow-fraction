#include "fileio.h"

#include<fstream>

#include "compression.h"
#include "date.h"
#include "mrtcalculator.h"
#include "struct.h"

#define SEGMENTED_BIN "seg"
#define DEPTH_IMAGE "depth_filtered_compressed"
#define GSV_IMAGE "gsv"
/*
#define SEGMENTED_BIN_POSTFIX "_0_0_seg.bin"
#define SEGMENTED_BIN_POSTFIX_LEN 12
*/

#define SEGMENTED_BIN_POSTFIX "_0_0.bin"
#define SEGMENTED_BIN_POSTFIX_LEN 8

#define DEPTH_IMAGE_POSTFIX "depth.png"

#define SEGMENTED_BIN_SIZE (512*512)
//#define DEPTH_IMAGE_SIZE (1048576)
#define DEPTH_IMAGE_SIZE (512*512)

using std::string;
using file::decompress;
namespace fs = std::experimental::filesystem;

namespace file {
	FileIO::FileIO(const JsonStruct *json) {
		inputDir = json->iRoot;
		outputDir = fs::path(json->oRoot).string();
		std::cout << "output directory: " << outputDir << std::endl;
		this->json = json;
		mrt.init(json->timeStep);
		/*if (fisheye == nullptr) fisheye = new Fisheye((GLFWwindow*)mrt.scSuf->getWindow());
		else {
			delete fisheye;
			fisheye = new Fisheye((GLFWwindow*)mrt.scSuf->getWindow());
		}*/

		fisheye_buffer.resize(6, nullptr);
		for (int i = 0; i < 6; i++) {
			fisheye_buffer[i] = (unsigned char*)malloc(SEGMENTED_BIN_SIZE * sizeof(unsigned char));
		}

		fisheye = new Fisheye();

		output.clear();
		for (int i = 0; i < 6; i++) {
			//TODO: when using a binary file, the count of channels would be one since we only use the red channel for depth image.
			unsigned char* dataPointer = (unsigned char*)malloc(DEPTH_IMAGE_SIZE * sizeof(unsigned char));
			depth_data.push_back(dataPointer);
			dataPointer = (unsigned char*)malloc(SEGMENTED_BIN_SIZE * sizeof(unsigned char));
			segmented_data.push_back(dataPointer);
			viewFactor.push_back({ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 });
		}		
#ifdef _DEBUG
		fs::create_directories("./Debug_output");
		fs::create_directories("./Debug_output/fisheye");
		fs::create_directories("./Debug_output/nonshaded_fisheye");
		fs::create_directories("./Debug_output/shaded_fisheye");
		fs::create_directories("./Debug_output/accumulated_temperature");
		debug_image_buffer = (unsigned char*)malloc(SEGMENTED_BIN_SIZE * 4 * sizeof(unsigned char));
		ilInit();
#endif // DEBUG
	}

	void FileIO::load() {
		calculate::Date startDate(json->startTime);
		calculate::Date endDate(json->endTime);
		calculate::Date temp(startDate);
		calculate::MRTCalculator calcmrt;
		std::cout << "start processing data" << std::endl;
		int gap_minute = endDate - startDate;
		unsigned len = (endDate - startDate)/json->timeStep + 1;
		std::cout << "total time steps: " << len << std::endl << std::endl;
		if (len != json->altitude->size() 
			|| len != json->airTemperature->size() 
			|| len != json->relativeHumidity->size() 
			|| len != json->windVelocity->size() 
			|| len != json->cloudCover->size()) {
			std::cout << "Incooperate length of meteorological data with time steps!\
 time steps: " << len << std::endl;
			exit(1);
		}
		for (int i = 0; i < len; i++) {
			temp.getDayOfYear();
			calculate::Date *t = new calculate::Date(temp);
			t->to_string();
			std::string *s = new string;
			output.push_back({ t, s });
			temp += json->timeStep;
		}

		//readingDir = inputDir;
		//genDir(readingDir, 0);
		std::string baseDir = "";
		std::string tileDir = "";
		recGenDir(inputDir, baseDir, tileDir, false);
	}

	void FileIO::write() {
		std::ofstream myfile;

		for (int i = 0; i < output.size(); i++) {
			std::string outname = outputDir + "\\" + output[i].first->str + writingDir;
			fs::create_directories(outname);
			outname += "\\result.csv";
			std::cout << outname << std::endl;
			myfile.open(outname, std::ofstream::out);
			myfile << csv_col;
			myfile << *(output[i].second);
			*(output[i].second) = "";

			myfile.close();
		}
	}

#ifdef _DEBUG
	void FileIO::saveImage(const unsigned char *data, std::string location, int size) {
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
#endif //_DEBUG

	void FileIO::recGenDir(const std::string& instring, const std::string& baseDir, const std::string& tileDir, bool flag) {
		fs::path inpath(instring);
		std::string name = inpath.filename().string();
		std::string base = baseDir;
		std::string tile = tileDir;
		bool state = flag;

		if (name == SEGMENTED_BIN) {
			state = true;
		}
		else if (name == DEPTH_IMAGE) {
			return;
		}
		else if (name == GSV_IMAGE || name == "seg_color" || name == "depth_filtered") {
			return;
		}
		else {
			if (state) {
				tile += "/" + name;
			}
			else {
				if (base.size() == 0) {
					base = inputDir;
				}
				else {
					base += "/" + name;
				}
			}
		}

		int res = 0;
		for (auto &dir : fs::directory_iterator(inpath)) {
			std::string input = dir.path().string();
			if (input.find(SEGMENTED_BIN_POSTFIX) != std::string::npos) {
				writingDir = tile;
				processData(input, base, tile);
				res++;
			} else recGenDir(input, base, tile, state);
		}
		
		if (res > 0) {
			write();
		}
	}

	int FileIO::loadDepth(std::vector<unsigned char*>& data, const std::string& depthDir) {

		/*for (int i = 0; i < 6; i++) {
			std::string dir = depthDir + depth_postfix[i];
			if (IL_TRUE == ilLoadImage(&dir[0])) {
				unsigned char* temp = ilGetData();
				for (int j = 0; j < DEPTH_IMAGE_SIZE; j++)
					data[i][j] = temp[j];
			}
			else {
				std::cout << "error: failed to load file: " << dir << std::endl;
			}
			
			//ILint size = ilGetInteger(IL_IMAGE_SIZE_OF_DATA);
		}*/
		for (int i = 0; i < 6; i++) {
			std::string dir = depthDir + depth_postfix[i];
			if (decompress(dir, data[i], SEGMENTED_BIN_SIZE) != 0) {
				std::cout << "failed when decompressing " << dir << std::endl;
				return 1;
			}
		}
		return 0;
	}

	int FileIO::loadSegmentedBin(std::vector<unsigned char*>& data, const std::string& segmentedDir) {
		for (int i = 0; i < 6; i++) {
			std::string dir = segmentedDir + fe_postfix[i];
			if (decompress(dir, data[i], SEGMENTED_BIN_SIZE) != 0) {
				std::cout << "failed when decompressing " << dir << std::endl;
				return 1;
			}
		}
		return 0;
	}

	void FileIO::processData(const std::string& instring, const std::string& baseDir, const std::string& tileDir) {
		fs::path inpath(instring);
		std::string name = inpath.filename().string();
		std::string latLon = name.substr(0, name.length() - SEGMENTED_BIN_POSTFIX_LEN);

		std::string depthDir = baseDir + "/" + DEPTH_IMAGE + tileDir + "/" + latLon;
		if (loadDepth(depth_data, depthDir)) return;

		std::string segDir = baseDir + "/" + SEGMENTED_BIN + tileDir + "/" + latLon;
		if (loadSegmentedBin(segmented_data, segDir)) return;
		fisheye->getVF(segDir, fisheye_buffer, vfCalculator, viewFactor);

#ifdef _DEBUG
		for (int i = 0; i < 6; i++) {
			//file::compress(fisheye_buffer[i], "./Debug_output/fisheye/" + latLon + postfix[i], 512 * 512, Z_DEFAULT_COMPRESSION);
			saveImage(fisheye_buffer[i], "./Debug_output/fisheye/" + latLon + debug_postfix[i], 512 * 512);
		}
#endif //DEBUG

		std::string tile = inpath.string().substr(
			inpath.parent_path().parent_path().string().length() + 1,
			inpath.parent_path().string().length() - inpath.parent_path().parent_path().string().length() - 1
		);

		GeoStruct geo;

		if (latLon.find('_') == string::npos) {
			std::cerr << "can't identify the LatLng of: " << name << std::endl;
		}
		else {
			int index = latLon.find('_');
			geo.lat = stod(latLon.substr(0, index));
			geo.lng = stod(latLon.substr(index + 1, latLon.length() - 1 - index));
		}

		if (tile.find('_') == string::npos) {
			std::cerr << "can't identify tile of: " << tile << std::endl;
		}
		else {
			tile[tile.find('_')] = ',';
		}

		mrt.pre_day = -1;

		for (int j = 0; j < output.size(); j++) {
			geo.altitude = json->altitude->at(j);
			geo.airTemperature = json->airTemperature->at(j);
			geo.windVelocity = json->windVelocity->at(j);
			geo.cloudCover = json->cloudCover->at(j);
			geo.relativeHumidity = json->relativeHumidity->at(j);

			*(output[j].second) += tile + ",";
			mrt.calculate(*(output[j].first), fisheye_buffer, geo, viewFactor[0], output[j].second, viewFactor, json, depth_data, fisheye, vfCalculator, segmented_data);
		}
	}

	/*void FileIO::genDir(std::string &inpath, int level) {
		fs::path filepath(inpath);
		if (level == 0) {
			for (auto &dir : fs::directory_iterator(filepath)) {
				if (fs::is_directory(dir)) {
					writingDir = outputDir + "\\" + dir.path().filename().string();
					for (const auto &i : output) {
						std::string dir = writingDir + "\\" + i.first->str;
						//std::cout << "creating folder: " << dir << std::endl;
						fs::create_directories(dir);
					}
					genDir(dir.path().string(), 1);
				}
			}
		}
		else if (level == 2) {
			for (auto &dir : fs::directory_iterator(filepath)) {
				if (fs::is_directory(dir)) {
					tileGroup = dir.path().filename().string();
					std::ofstream myfile;
					for (int i = 0; i < output.size(); i++) {
						std::string outname = writingDir + "\\" + output[i].first->str + "\\" + tileGroup;
						std::cout << "output folder name: " << outname << std::endl;
						fs::create_directory(outname);
						outname += "\\result.csv";
						myfile.open(outname);

						myfile << csv_col;

						myfile.close();
					}

					recursiveLoad(dir.path().string(), true);

					for (int i = 0; i < output.size(); i++) {
						std::string outname = writingDir + "\\" + output[i].first->str + "\\" + tileGroup;
						outname += "\\result.csv";
						myfile.open(outname, std::ofstream::out | std::ofstream::app);
						myfile << std::endl;

						myfile.close();
					}
				}
			}
		}
		else if (level == 1) {
			for (auto &dir : fs::directory_iterator(filepath)) {
				if (fs::is_directory(dir)) {
					genDir(dir.path().string(), 2);
				}
			}
		}
	}

	void FileIO::write() {
		std::ofstream myfile;
		
		for (int i = 0; i < output.size(); i++) {
			std::string outname = writingDir + "\\" + output[i].first->str + "\\" + tileGroup;
			outname+= "\\result.csv";
			myfile.open(outname, std::ofstream::out | std::ofstream::app);

			myfile << *(output[i].second);
			*(output[i].second) = "";

			myfile.close();
		}
	}*/

	/* recursively read directories and read image or binary files
	to generate fisheye file. Para basepath is the base path of input.
	Inpath is the path of the current directory. Output is the output string.
	Deprecated: 10/9/2018
	*/
	int FileIO::recursiveLoad(const string& inpath, bool group) {
		fs::path filepath(inpath);
		string inDir;
		fs::path base(inputDir);
		inDir = base.string();

		string curDir = filepath.string();

		int res = 0;
		int total_dir = 0;
		for (auto &dir : fs::directory_iterator(filepath)) {
			if (fs::is_directory(dir)) {
				total_dir++;
			}
		}

		int dir_count = 0;

		for (auto &dir : fs::directory_iterator(filepath)) {
			if (fs::is_directory(dir)) {
				string temp = dir.path().string();
				std::clock_t start;
				double duration;
				start = std::clock();
				int count = recursiveLoad(temp, false);
				duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
				dir_count++;
				std::cout << "Finished processing files in folder " << temp << std::endl;
				std::cout << ". Time used: " << duration
					<< ". Files processed: " << count << ". "
					<< "Folder processed: " << dir_count << "/" << total_dir
					<< std::endl;
				continue;
			}

			string s = dir.path().string();
			if (s.find("_0_90.bin") != string::npos) {
				string file_name_path = s.substr(
					0, s.length() - bin_postfix[s.substr(s.length() - 8, 8)]
				);

				string name = s.substr(
					dir.path().parent_path().string().length(),
					s.length() - 13 - dir.path().parent_path().string().length()
				);

				string tile = s.substr(
					filepath.parent_path().string().length()+1,
					dir.path().parent_path().string().length() - filepath.parent_path().string().length()-1
				);

				unsigned char *data = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));
				GeoStruct geo;
				
				if (name.find('_') == string::npos) {
					std::cerr << "can't identify the LatLng of: " << name << std::endl;
				}
				else {
					int index = name.find('_');
					geo.lat = stod(name.substr(1, index - 1));
					geo.lng = stod(name.substr(index+1, name.length() - 1 - index));
				}
				
				if (tile.find('_') == string::npos) {
					std::cerr << "can't identify tile of: " << tile << std::endl;
				}
				else {
					tile[tile.find('_')] = ',';
				}

				/*for (int i = 0; i < 6; i++) {
					if (i == 2) continue;
					string data_path = file_name_path + postfix[i];
					fisheye->
					decompress(data_path, data, 512 * 512);
					vfCalculator.calculate(data, viewFactor[i]);
				}

				string data_path = file_name_path + postfix[2];
				decompress(data_path, data, 512 * 512);
				vfCalculator.calculate(data, viewFactor[2]);
				*/
				/*fisheye->getVF(file_name_path, data, vfCalculator, viewFactor);

				for (int j = 0; j < output.size(); j++) {
					geo.altitude = json->altitude->at(j);
					geo.airTemperature = json->airTemperature->at(j);
					geo.windVelocity = json->windVelocity->at(j);
					geo.cloudCover = json->cloudCover->at(j);
					geo.relativeHumidity = json->relativeHumidity->at(j);

					*(output[j].second) += tile + ",";
					mrt.calculate(*(output[j].first), data, geo, viewFactor[0], output[j].second, viewFactor, json);
				}
				delete[] data;
				res++;*/
			}
		}
		if (res > 0) write();

		return res;
	}

	void FileIO::dellocateBuffer(vector<unsigned char*> &buffer) {
		for (unsigned char *p : buffer) {
			if (p) free(p);
		}
		buffer.clear();
	}

	FileIO::~FileIO() {
		if (fisheye) delete fisheye;
		fisheye = nullptr;

		dellocateBuffer(depth_data);
		dellocateBuffer(segmented_data);
		dellocateBuffer(fisheye_buffer);

		for (auto &p : output) {
			delete p.first;
			delete p.second;
		}
		output.clear();
	}
} // namespace file