#ifndef CALCULATEMRT_STRUCT_H_
#define CALCULATEMRT_STRUCT_H_

namespace file {
	struct GeoStruct {
		double lat;
		double lng;
		
		double altitude;
		double cloudCover;
		double relativeHumidity;
		double airTemperature;
		double windVelocity;
	};
};


#endif // !CALCULATEMRT_STRUCT_H_