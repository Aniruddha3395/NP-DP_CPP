# NON-PLANAR 3D PRINTING TOOLPATH GENERATION

### lIBRARIES REQUIRED:
1. Eigen ([Link](http://eigen.tuxfamily.org/index.php?title=Main_Page))
2. Libigl ([Link](https://github.com/libigl/libigl)) for reading STL file
3. gnuplot-iostream ([Link](https://github.com/dstahlke/gnuplot-iostream)) for visualization data (optional)

### INSTALLATION INSTRUCTIONS FOR LIBIGL AND DEPENDANCIES:

1. download zip or clone git directory (https://github.com/libigl/libigl)
2. install all the dependancies listed 

```
 sudo apt-get install git
 sudo apt-get install build-essential
 sudo apt-get install cmake
 sudo apt-get install libx11-dev
 sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev
 sudo apt-get install libxrandr-dev
 sudo apt-get install libxi-dev
 sudo apt-get install libxmu-dev
 sudo apt-get install libblas-dev
 sudo apt-get install libxinerama-dev
 sudo apt-get install libxcursor-dev
```

2. build steps :-	mkdir build
			cd build/
			cmake ..
			make -j 
			sudo make install
3. NOTE : If only header files are added in the /usr/local/include/igl folder, then add all the cpp and extra support files from cloned libigl source directory (copy all files from ${SOURCE_DIR}/include/igl/ )
4. To avoid issue of 'pythread', add [target_link_libraries(${EXECUTABLE_FILE} pthread) to cmake]

### INSTALLATION INSTRUCTIONS FOR GNUPLOT-IOSTREAM AND DEPENDANCIES:

1. download zip or clone git directory ([Link](https://github.com/dstahlke/gnuplot-iostream))
2. include the directory in CMakeLists and header in main_file.cpp
3. certain dependancies required:
```
[target_link_libraries(${EXECUTABLE_FILE} ${PATH}/libboost_iostreams.a)]
[target_link_libraries(${EXECUTABLE_FILE} ${PATH}/libboost_system.a)]
[target_link_libraries(${EXECUTABLE_FILE} ${PATH}/libboost_filesystem.a)]
[target_link_libraries(${EXECUTABLE_FILE} ${PATH}/libutil.a)]
```
### RUNNING MAIN FILE:

1. add stl file in the data/
2. specify appropriate parameters in main_file.cpp
3. to generate the .out file-
   in the main directory-
```
mkdir build
cd build/
cmake ..
make -j
./OUTFile
```

 
