cd ../calculateMRT
mkdir -p ./build
cd ./build
cmake ..
cmake --build . --config Release
cd ..
cd ..
cd ./release
cp ../calculateMRT/build/src/Release/calculateMRT.exe ./
cp ../calculateMRT/bin/* ./
cp ../calculateMRT/config.json ./
cp ../calculateMRT/README.md ./
mkdir -p ./shaders
cp ../calculateMRT/shaders/* ./shaders/
rm -rf ../calculateMRT/build
