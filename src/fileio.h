#ifndef CALCULATEMRT_FILEIO_H_
#define CALCULATEMRT_FILEIO_H_

#include<experimental/filesystem>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<unordered_map>
#include<ctime>
#include<utility>
#include<vector>

#include"IL/il.h"

#include "jsonstruct.h"
#include "date.h"
#include "vfcalculator.h"
#include "mrtcalculator.h"
#include "fisheye/fisheye.h"

namespace file {
	/**
		@class FileIO
		@brief

		This is the entry of processing data. Loading segmented and depth files from folder structures. 
		Write calculate results to output directories.
	*/
	class FileIO {
	public:
		/**
			@brief Constructs the object and initialize input and output directory from json structure. 
			@param json Json struct including configuration
			Initialize mrt calculator and data buffers.
		*/
		FileIO(const JsonStruct *json);
		
		/**
			@brief Create time stamps and init output directory
		*/
		void load();

		/**
			@brief Read depth image files into buffer
			@param data Depth data buffer
			@param depthDir Directory of depth images
		*/
		int loadDepth(std::vector<unsigned char*> &data, const std::string& depthDir);

		/**
			@brief Read segmented binary files into buffer
			@param data Segmented binary data buffer
			@param segmentedDir Directory of segmented binary files
		*/
		int loadSegmentedBin(std::vector<unsigned char*> &data, const std::string& segmentedDir);
		/**
			@brief Output result.csv
		*/
		void write();

		// deprecated 10-9
		int recursiveLoad(const std::string& inpath, bool group);
		//void genDir(std::string &inpath, int level);

		/**
			@brief Recursively read directory parse base path and relative path
			@param inpath Current input directory
			@param baseDir Base directory which contains depth folder, segmented folder and image folder.
			@param tileDir Relative directory which contains tile directory
			@param flag Symbol to record current directory should belong to base or tile.
		*/
		void recGenDir(const std::string& inpath, const std::string& baseDir, const std::string& tileDir, bool flag);

		/**
			@brief Calculate result for one position
			@param inpath Input directory of data
			@param basDir Base directory of diverse data folders
			@param tileDir Relative directory of data folders
		*/
		void processData(const std::string& inpath, const std::string& baseDir, const std::string& tileDir);
		Fisheye *fisheye = nullptr;
		~FileIO();

	private:

		std::string inputDir;
		std::string outputDir;

		std::string readingDir;
		std::string writingDir;

		std::string tileGroup;
		const JsonStruct *json;

		//TODO: needs to deallocate memory
		std::vector<unsigned char*> depth_data;
		std::vector<unsigned char*> segmented_data;

		std::vector<unsigned char*> fisheye_buffer;
	
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

		/**
		@brief Save images for debugging purpose
		@param data Unsigned char pointer of image input
		@param location Directory to output file
		@param size of image data
		*/
		void saveImage(const unsigned char *data, std::string location, int size);
#endif //_DEBUG

/*		const std::string depth_postfix[6] = {
			"_0_0_depth.png",    //north
			"_180_0_depth.png",  //south
			"_270_0_depth.png", //east
			"_90_0_depth.png", //west
			"_0_90_depth.png",   //up
			"_0_270_depth.png"  //down
		};

		std::unordered_map<std::string, int> bin_postfix = {
			{ "_0_0.bin", 8 },
			{ "0_90.bin", 9 },
			{ "_270.bin", 10 },
			{ "90_0.bin", 9 },
			{ "80_0.bin", 10 },
			{ "70_0.bin", 10 }
		};

		std::string postfix[6] = {
			"_0_90_seg.bin",   //up
			"_0_270_seg.bin",  //down
			"_0_0_seg.bin",    //north
			"_180_0_seg.bin",  //south
			"_90_0_seg.bin",   //east
			"_270_0_seg.bin"  //west
		};

		const string fe_postfix[6] = {
			"_90_0_seg.bin",
			"_270_0_seg.bin",
			"_0_90_seg.bin",
			"_0_270_seg.bin",
			"_0_0_seg.bin",
			"_180_0_seg.bin"
		};*/
		const std::string depth_postfix[6] = {
			"_0_0.bin",    //north
			"_180_0.bin",  //south
			"_270_0.bin", //east
			"_90_0.bin", //west
			"_0_90.bin",   //up
			"_0_270.bin"  //down
		};

		std::unordered_map<std::string, int> bin_postfix = {
			{ "_0_0.bin", 8 },
			{ "0_90.bin", 9 },
			{ "_270.bin", 10 },
			{ "90_0.bin", 9 },
			{ "80_0.bin", 10 },
			{ "70_0.bin", 10 }
		};

		std::string postfix[6] = {
			"_0_90.bin",   //up
			"_0_270.bin",  //down
			"_0_0.bin",    //north
			"_180_0.bin",  //south
			"_90_0.bin",   //east
			"_270_0.bin"  //west
		};

		const string fe_postfix[6] = {
			"_90_0.bin",
			"_270_0.bin",
			"_0_90.bin",
			"_0_270.bin",
			"_0_0.bin",
			"_180_0.bin"
		};

    std::clock_t timer;

		// output buffer
		std::vector<std::pair<calculate::Date *, std::string *>> output;
		calculate::VFCalculator vfCalculator;
		//std::vector<double> viewFactor;
		calculate::MRTCalculator mrt;

		// view factor data buffer
		std::vector<std::vector<double>> viewFactor;

		// free memory
		void dellocateBuffer(vector<unsigned char*> &buffer);

		// column names of output .csv files
		const std::string csv_col = "Tile-x, Tile-y, Latitude, Longitude, Date, Day of Year,\
	Sun Hours, Sun Visible,\
	Sky_E, Tree_E_Exposed, Building_E_Exposed, Impervious_E_Exposed, Pervious_E_Exposed, Moving_E_Exposed, Tree_E_Shaded, Building_E_Shaded, Impervious_E_Shaded, Pervious_E_Shaded, Moving_E_Shaded,\
  Sky_W, Tree_W_Exposed, Building_W_Exposed, Impervious_W_Exposed, Pervious_W_Exposed, Moving_W_Exposed, Tree_W_Shaded, Building_W_Shaded, Impervious_W_Shaded, Pervious_W_Shaded, Moving_W_Shaded,\
  Sky_Up, Tree_Up_Exposed, Building_Up_Exposed, Impervious_Up_Exposed, Pervious_Up_Exposed, Moving_Up_Exposed, Tree_Up_Shaded, Building_Up_Shaded, Impervious_Up_Shaded, Pervious_Up_Shaded, Moving_Up_Shaded,\
  Sky_Down, Tree_Down_Exposed, Building_Down_Exposed, Impervious_Down_Exposed, Pervious_Down_Exposed, Moving_Down_Exposed, Tree_Down_Shaded, Building_Down_Shaded, Impervious_Down_Shaded, Pervious_Down_Shaded, Moving_Down_Shaded,\
  Sky_N, Tree_N_Exposed, Building_N_Exposed, Impervious_N_Exposed, Pervious_N_Exposed, Moving_N_Exposed, Tree_N_Shaded, Building_N_Shaded, Impervious_N_Shaded, Pervious_N_Shaded, Moving_N_Shaded,\
  Sky_S, Tree_S_Exposed, Building_S_Exposed, Impervious_S_Exposed, Pervious_S_Exposed, Moving_S_Exposed, Tree_S_Shaded, Building_S_Shaded, Impervious_S_Shaded, Pervious_S_Shaded, Moving_S_Shaded,\
  Sky_Total, Shadow_Total, Exposed_Total\
  \n";
	};
}

#endif // !CALCULATEMRT_FILEIO_H_