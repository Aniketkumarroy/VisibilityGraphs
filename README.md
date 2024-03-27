# VisibilityGraphs

## Build
for linking and building the library with your source code add your file in [CMakelists.txt](CMakeLists.txt) in line no. 9

then as usual
create a build directory and navigate into it
```bash
mkdir build
cd build
```
now invoke cmake and make commands
```bash
cmake ..
make
```
your executable will be ready in build directory with the same name as your source file name without any extension

## Use
for getting all the visible points from point ```p``` by using the visibility graph created using ```vg.build()``` use 

    std::vector<Point> visible = vg.VG_Ptr->get_adjacent_points(p);
it returns a vector of all the visible points.

if you have not build the visibility graph but want visible vertices for a point ```p``` and graph ```G``` use 

    std::unordered_set<Point> visible = vg.find_visible_vertices(p, G);
it returns a std::unordered_set of all the visible points.