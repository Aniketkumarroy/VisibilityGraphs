#ifndef VIS_GRAPH
#define VIS_GRAPH

#include"visible_vertices.h"
#include"ThreadPool.h"

class VisGraph{
    public:
        VisGraph();

        void SingleThreadProcessing(const std::vector<Point>& points, int& batch_size);

        void MultiThreadProcessing(const std::vector<Point>& points, int& batch_size, int& threads);

        void build(std::vector<std::vector<Point>> &polygons, std::pair<int,int> grid_resolution = {-1, -1}, int batch_size=10, int threads=1, bool status=false);

        void UpdateByADD(std::vector<std::vector<Point>> &polygons, int& batch_size, int& threads, bool tolerance=false);

        std::unordered_set<Point> find_visible_vertices(const Point& point, const Graph& graph);

        Graph& graph();
        Graph& visgraph();
        void clear();


    public:
        std::unique_ptr<Graph> G_Ptr, VG_Ptr;
        bool graph_built;
};

#endif