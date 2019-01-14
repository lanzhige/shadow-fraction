file(GLOB hdrs 
	"compression.h"
	"jsonstruct.h" 
	"jsonreader.h"
	"date.h"
	"fileio.h"
	"mrtcalculator.h"
	"vfcalculator.h"
	"struct.h"
	"./fisheye/*.h"
  "fractioncalculator.h"
)

file(GLOB srcs 
	"fraction.cpp"
	"jsonreader.cc"
	"jsonstruct.cc"
	"date.cc"
	"fileio.cc"
	"mrtcalculator.cc"
	"vfcalculator.cc"
  "fractioncalculator.cc"
)

add_executable(shadow-fraction ${hdrs} ${srcs})

target_link_libraries(shadow-fraction ${Boost_LIBRARIES} OpenGL32.lib glew32s.lib glfw3.lib DevIL.lib zlib.lib debug SfcCalcd.lib optimized SfcCalc.lib)

install(TARGETS shadow-fraction DESTINATION .)