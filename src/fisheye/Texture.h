#pragma once

#include "../compression.h"
// GL Includes
#define GLEW_STATIC

#define IMG_SIZE 262144

#include <vector>
#include <string>
#include <fstream>

using std::vector;
using std::string;
using std::ifstream;

class TextureLoading {
public:    
	/*static GLuint loadBinary(const vector<string> &faces) { 
		GLuint textureID;
		glGenTextures(1, &textureID);

		int imageWidth = 512, imageHeight = 512;

		glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);
		char *tempBuffer = new char[IMG_SIZE];
		for (GLuint i = 0; i < faces.size(); i++) {
			ifstream myFile(faces[i], std::ios::in | std::ios::binary);
			if (!myFile.read(tempBuffer, IMG_SIZE)) {
				std::cerr << "Reading file failed: " << faces[i] << std::endl;
				std::cerr << "read file status: " << myFile.gcount();
			}
			myFile.close();
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_R8, imageWidth, imageHeight, 0, GL_RED, GL_UNSIGNED_BYTE, tempBuffer);
		}
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
		delete[] tempBuffer;
		return textureID;
	}*/
	static GLuint loadBinary(const vector<string> &faces) {
		GLuint textureID;
		glGenTextures(1, &textureID);

		int imageWidth = 512, imageHeight = 512;

		glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);
		//char *tempBuffer = new char[IMG_SIZE];
		unsigned char *data = (unsigned char *)malloc(512 * 512 * sizeof(unsigned char));

		for (GLuint i = 0; i < faces.size(); i++) {
			/*ifstream myFile(faces[i], std::ios::in | std::ios::binary);
			if (!myFile.read(tempBuffer, IMG_SIZE)) {
			std::cerr << "Reading file failed: " << faces[i] << std::endl;
			std::cerr << "read file status: " << myFile.gcount();
			}*/
			file::decompress(faces[i], data, 512 * 512);
			//myFile.close();

			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_R8, imageWidth, imageHeight, 0, GL_RED, GL_UNSIGNED_BYTE, data);
		}
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
		//delete[] tempBuffer;
		delete[] data;
		return textureID;
	}

	static GLuint loadBinary(vector<unsigned char *>& segmented) {
		GLuint textureID;
		glGenTextures(1, &textureID);

		int imageWidth = 512, imageHeight = 512;

		glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

		for (GLuint i = 0; i < segmented.size(); i++) {
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_R8, imageWidth, imageHeight, 0, GL_RED, GL_UNSIGNED_BYTE, segmented[i]);
		}
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_CUBE_MAP, 0);
		//delete[] tempBuffer;
		return textureID;
	}
};
