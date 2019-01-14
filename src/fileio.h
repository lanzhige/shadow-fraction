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
		const std::string csv_col = "MapTile-x, MapTile-y, Latitude, Longitude, Year, Month, Day, Hours, Day of Year, Air Temperature,\
	Relative Humidity,	Vapor Pressure, Wind Speed,	Sun Hours, Sun Visible,\
	Gamma, Phi, Pressure at given Altitude(Atmospheric), Relative Optical Air Mass,\
	Optical Thickness, Extra-Terrestrial Radiation, Direct Radiation, Global Radiation,\
	Diffused Radiation,	Atmospheric Radiation, \
	SurfaceTemperatureSky, SurfaceTemperatureTree, SurfaceTemperatureBuilding, SurfaceTemperatureImpervious, SurfaceTemperaturePervious, SurfaceTemperatureMovingObjects, LongWaveRadiationSky, LongWaveRadiationTree, LongWaveRadiationBuilding, LongWaveRadiationImpervious, LongWaveRadiationPervious, LongWaveRadiationMovingObjects, \
	SVF_Up, TVF_Up,	BVF_Up,	IVF_Up, PVF_Up, MOVF_Up, \
	SVF_Down, TVF_Down,	BVF_Down, IVF_Down,	PVF_Down, MOVF_Down, \
	SVF_N, TVF_N, BVF_N, IVF_N,	PVF_N, MOVF_N, \
	SVF_S, TVF_S, BVF_S, IVF_S, PVF_S, MOVF_S, \
	SVF_E, TVF_E, BVF_E, IVF_E, PVF_E, MOVF_E, \
	SVF_W, TVF_W, BVF_W, IVF_W,	PVF_W, MOVF_W, \
	SVF_360, TVF_360, BVF_360, IVF_360, PVF_360, MOVF_360, SVF_MRT, TVF_MRT, IVF_MRT, BVF_MRT, PVF_MRT, MOVF_MRT, \
    MRT_North, MRT_Down, MRT_Up, MRT_South, MRT_West, MRT_East, \
	Mean Radiant Temperature(MRT), \
	UTCI, PET, Core Temperature(C), Skin Temperature(C), Total Water Loss(g / h), Skin Wetting, Internal Heat(W), Radiation Balance(W), Convection(W), Water Vapor Diffusion(W), Welding Evaporation(W), Respiration(W), \
	SD_TVF_Up, SD_BVF_Up,SD_IVF_Up,SD_PVF_Up,SD_MOVF_Up, \
	SD_TVF_Down, SD_BVF_Down, SD_IVF_Down,	SD_PVF_Down, SD_MOVF_Down, \
	SD_TVF_N, SD_BVF_N, SD_IVF_N, SD_PVF_N, SD_MOVF_N, \
	SD_TVF_S, SD_BVF_S, SD_IVF_S, SD_PVF_S, SD_MOVF_S, \
	SD_TVF_E, SD_BVF_E, SD_IVF_E, SD_PVF_E, SD_MOVF_E, \
	SD_TVF_W, SD_BVF_W, SD_IVF_W, SD_PVF_W, SD_MOVF_W, \
	NSD_TVF_Up, NSD_BVF_Up, NSD_IVF_Up, NSD_PVF_Up, NSD_MOVF_Up, \
	NSD_TVF_Down, NSD_BVF_Down, NSD_IVF_Down, NSD_PVF_Down, NSD_MOVF_Down, \
	NSD_TVF_N, NSD_BVF_N, NSD_IVF_N, NSD_PVF_N, NSD_MOVF_N, \
	NSD_TVF_S, NSD_BVF_S, NSD_IVF_S, NSD_PVF_S, NSD_MOVF_S, \
	NSD_TVF_E, NSD_BVF_E, NSD_IVF_E, NSD_PVF_E, NSD_MOVF_E, \
	NSD_TVF_W, NSD_BVF_W, NSD_IVF_W, NSD_PVF_W, NSD_MOVF_W\n";
	};
}

#endif // !CALCULATEMRT_FILEIO_H_