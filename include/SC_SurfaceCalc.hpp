#pragma once

/*! 

	@mainpage SurfaceCalc - A Library for the caculation of surface features

	@section Overview

	SurfaceCalc is a library, which calculates normals and shade information based on cubic depth images and the sun position.
	It projects the depth images back into 3D space and reconstructs the geometry by triangularization. The resulting mesh
	is then used to extract the surface normals. With mesh, normals and the position of the sun, the shade information can be calculated 
	through shadow mapping.

	SurfaceCalc is implemented in C++ and uses the OpenGL library to achieve a rough performance of 0.1 sec/location, where each
	one contains 6 depth images. It was developed on and for Windows 10.

	Additional parameters allow the debugging and exploration of the projected depth images in realtime.

	The calculation is done by projecting the depth images into 3D space while triangulizing the pixels, resulting in an deformed plane
	mesh. Since the distortion can be very high, some triangles need to be discarded. The "constant DDT" and "projective DDT" are thresholds
	(Depth Difference Threshold) upon which a triangle is discarded. If the difference in depth to one of its neighbors exceeds the threshold
	the triangle is not rendered, as it is considered to be an artifact of the projection and not actually part of the reconstructed environment.
	The differentiation between constant and projective DDT comes from the way depth is handled in modern rendering pipelines. Usually
	the depth buffer stores the values non linearly, so that greater precision is provided for nearer objects, which when ordered wrong would
	result in an much greater an obvious error than distant objects. Therefore the precision of the depth buffer is often adjusted like depth=1/z.
	It is therefore important to provide an additional threshold for this phenomenon. The projective DDT is then scaled with the w component of
	the projected vertex to mimic the distribution of the depth buffer. Unfortunately the constant and projective threshold cannot be combined
	since the actual depth function is unknown. The combination of the threshold has, however, empirically proven to be effective.

	@subsection Depthimages

	ShadowCalc uses six depth images, which are aligned like a cube to reconstruct the environment. These depth images have therefore
	a FOV of 90 degrees and a aspect ratio of 1. Although these assumptions are not required by the tool, they produce the best results,
	since they then would map the hole environment.

	The depth images have to be form the outside, so that the sun can be simulated an the shade information is meaningful.
	To keyout the sky, the depth of all skypixels has to be set to the maximum (255 for unsigned char).

	@subsection Surfaceimages

	The resulting surface images are returned in RGBA float row-major format and contain the information of (normal.x, normal.y, normal.z, shade),
	where the shade information is 1 if the pixel is in the shade and 0 if not.

	@section Dependencies

	The main library utilizes the following libraries:

	- Glew library (http://glew.sourceforge.net/)

	- GLFW library (https://www.glfw.org/)

	- OpenGL (https://www.opengl.org/)

	- GLM library (https://glm.g-truc.net/0.9.9/index.html)

	Since GLM is a header only library, it is not needed in order to use SurfaceCalc. The others must be linked to, though.
	
	The exmaple uses an additionally library for loading and storing images:

	- DevIL image loading library (http://openil.sourceforge.net/)

	 @section Guide

	 The usage of the library is straight forward. Add the /include folder to your include paths and link to the 
	 static library found in /bin. Additionally you have to link to GLEW, GLFW and OpenGL libraries has well. 
	 A copy of each one is found in the /lib folder.

	 A working example is found in the /example folder, which also illustrates the appropriate usage of the tool.

	 @section References

	 - The depth images found in /examples/data were estimated based on color images by "monodepth" (https://github.com/mrharicot/monodepth)

	 - The color images found in /examples/data are from Google Street View and were provided by Prof. Dr. Ariane Middel  

	 - The segmentation images found in /examples/data were derived form another project in association with Prof. Dr. Ariane Middel

	 @section License

	 << in development >>
	 
*/

//other includes
#include <vector>

//forward declaration
class SC_SfcApp;

/**
	@class SC_SurfaceCalc
	@brief

	Is the central interface for the general use of this tool. All methods used form outside go through here. All parameters, that define the 
	behaviour are also located here and can easily be accessed or changed. There should only be one interface per program. The behaviour of multiple
	instances of this class is UNDEFINED.
*/
class SC_SurfaceCalc{

public:

	/**
		@brief Constructs the interface and intializes all the memory associated with it. 
	*/
	SC_SurfaceCalc();

	/**
		@brief Destroys the interface and frees all the memory associated with it.
	*/
	~SC_SurfaceCalc();

	/**
		@brief Initializes the state of the program to match the adjusted parameters. !! Has to be called before other methods can be used !!
	*/
	void initialize();

	/**
		@struct Parameters
		@brief

		Is a combination of all variables, taht define the behaviour of the program. They can be changed from outside, but are initialized to
		meaningful values. The changes take only effect before initialize is called.
	*/
	struct Parameters{

		//pass optionally data to display color and segmentation when in window mode
		std::vector<unsigned char*> DATA_SEGMENTATION;
		std::vector<unsigned char*> DATA_COLOR;

		//shader options
		float INTERPOLATION_FACTOR_SEGMENTATION = 0.5f;
		float INTERPOLATION_FACTOR_COLOR		= 0.5f;
		float INTERPOLATION_FACTOR_DIFFUSE		= 0.5f;
		float INTERPOLATION_FACTOR_SHADOW		= 0.5f;
		bool ENABLE_KEYOUT_SKY					= true;
		bool ENABLE_DISPLAY_NORMALS				= false;
		bool ENABLE_DISPLAY_SUN					= false;
		float DDT_CONSTANT_CULLING_FACTOR		= 0.01f;
		float DDT_PROJECTIVE_CULLING_FACTOR		= 0.02f;

		//shadowmap options
		int SHADOW_MAP_WIDTH	= 4096;
		int SHADOW_MAP_HEIGHT	= 4096;
		bool SHADOW_MAP_EXPORT	= false;

		//specify the properties of the depth images passed
		unsigned int INPUT_IMAGE_WIDTH				= 512;
		unsigned int INPUT_IMAGE_HEIGHT				= 512;
		unsigned int INPUT_IMAGE_COMPONENT_COUNT	= 4;

		//specify the projection of the vertices
		float VERTEX_INVERSE_PERSP_PROJ_FOV			= 90.0f;
		float VERTEX_INVERSE_PERSP_PROJ_ASPECT		= 1.0f;
		float VERTEX_INVERSE_PERSP_PROJ_NEAR		= 0.1f;
		float VERTEX_INVERSE_PERSP_PROJ_FAR			= 1000.0f;

		//specify the camera for the shadow mapping
		float SHADOW_ORTHO_PROJ_LEFT	= -7.0f;
		float SHADOW_ORTHO_PROJ_RIGHT	=  7.0f;
		float SHADOW_ORTHO_PROJ_TOP		= -7.0f;
		float SHADOW_ORTHO_PROJ_BOTTOM	=  7.0f;
		float SHADOW_ORTHO_PROJ_NEAR	= 1.0f;
		float SHADOW_ORTHO_PROJ_FAR		= 50.0f;

		//specify the camera in window mode
		float VIEW_PERSP_PROJ_FOV		= 60.0f;
		float VIEW_PERSP_PROJ_ASPECT	= 16.0f / 9.0f;
		float VIEW_PERSP_PROJ_NEAR		= 0.001f;
		float VIEW_PERSP_PROJ_FAR		= 5000.0f;
		float VIEW_MOVEMENT_SPEED		= 0.001f;
		float VIEW_ROTATION_SPEED		= 0.1f;

		//specify window porperties (need to be adjusted before initialization to take effect)
		bool DISPLAY_OUTPUT_IN_WINDOW		= false;
		unsigned int WINDOW_INITIAL_WIDTH	= 720;
		unsigned int WINDOW_INITIAL_HEIGHT	= 480;
		char *WINDOW_TITLE					= "SurfaceCalc";

	} parameters;

	/**
		@struct Screenshot
		@brief Is an output structure to export the shadow map or other screenshots of the application. Mainly used for debugging.
	*/
	struct Screenshot{

		int width;
		int height;
		unsigned char* data;

	};
	std::vector<Screenshot*> screenshots;

	/**
		@brief Function to comfortably delete all resources associated with a screenshot
		@param screenshot The screenshot, that is going to be deleted
	*/
	static void deallocateScreenshot(Screenshot *&screenshot);

	/**
		@brief Calcualtes the surface images of the given depth images and sun position
		@param depthImages Raw data of the depth images. Expects to be in the sequence of (north, south, east, west, up, down) and in row major format. It has also to comply with the parameters (e.g. INPUT_IMAGE_WIDTH)
		@param sunZenithAngle The angle between the sun and the horizon in degrees
		@param sunAzimuthAngle The angle between the sun and north in degrees
		@return Returns the raw data of the surface images in row major order and in sequence of (north, south, east, west, up, down). Each pixel contains the values (normal.x, normal.y, normal.z, shadow), value 0 indicating shadow and 1 indicating light (no values in between)
	*/
	std::vector<float*>* calculateSurfaceImages(std::vector<unsigned char*> &depthImages, const double sunZenithAngle, const double sunAzimuth) const;

	/**
		@brief Frees all memory associated with the surface images and clears the vector for later reuse.
		@param surfaceImages The list of surface images to be cleared
	*/
	static void deallocateSurfaceImages(std::vector<float*> *surfaceImages);

	/**
		@brief return the glfw window pointer
	*/
	void *getWindow();

private:

	//pointer to implementation (pimp)
	SC_SfcApp *m_application;

};