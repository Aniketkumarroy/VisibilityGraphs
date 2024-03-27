/*
    For compiling the code use the '-pyhread' flag as it uses std::thread for parallelism

    g++ -pthread main.cpp -o main

*/


#include<iostream>
#include<vector>
#include"visgraph.h"

int main(){

    Point p1(1028, 511.0);
    Point p2(853.0, 330.0);
    Point p3(1037.0, 204.0);
    Point p4(1185.0, 504.0);
    Point p5(124.0, 663.0);
    Point p6(68.0, 440.0);
    Point p7(292.0, 340.0);
    Point p8(544.0, 624.0);

    std::vector<std::vector<Point>> v;
    v.push_back({p1, p2, p3, p4});
    v.push_back({p5, p6, p7, p8});

    VisGraph vg;
    std::pair<int,int> grid_resolution = {1,1};
    int no_threads = 1, batch_size = 1;

    vg.build(v, grid_resolution, batch_size, no_threads);

    // for getting all the visible points from point 'p' by using the visibility graph breated using vg.build() 
    // use vg.VG_Ptr->get_adjacent_points(p), it returns a vector of all the visible points.
    // if you have not build the visibility graph but want visible vertices for a point 'p' and graph 'G' use
    // vg.find_visible_vertices(p, G), it will return an std::unordered_set<Point> of the visible points

    // print the graph
    std::cout<<vg.VG_Ptr->toString();
    return 0;
}