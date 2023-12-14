# VisibilityGraphs

For compiling the code use the ```-pyhread``` flag as it uses std::thread for parallelism

    g++ -pthread main.cpp -o main


for getting all the visible points from point ```p``` by using the visibility graph created using ```vg.build()``` use 

    std::vector<Point> visible = vg.VG_Ptr->get_adjacent_points(p);
it returns a vector of all the visible points.

if you have not build the visibility graph but want visible vertices for a point ```p``` and graph ```G``` use 

    std::unordered_set<Point> visible = vg.find_visible_vertices(p, G);
it returns a std::unordered_set of all the visible points.