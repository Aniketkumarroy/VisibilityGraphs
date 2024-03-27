#ifndef VIS_GRAPH
#define VIS_GRAPH

#include"visible_vertices.h"
#include"ThreadPool.h"

class VisGraph{
    public:
        VisGraph() : G_Ptr(nullptr), VG_Ptr(nullptr), graph_built(false) {}

        void SingleThreadProcessing(const std::vector<Point>& points, int& batch_size){
            int s = points.size();
            // progress = (status)? std::make_unique<ProgressBar>(20, s/batch_size): nullptr;
            for(int i=0; i<s; i=i+batch_size){
                int start = i;
                int end = (i + batch_size >= s)? s: i+batch_size;
                for(int j=start; j<end; j++){
                    for(auto p: visible_vertices(points[j], *G_Ptr, &points))
                    VG_Ptr->add_edge(Edge(p, points[j]));
                }
                // if (progress != nullptr) progress->update();
            }
        }

        void MultiThreadProcessing(const std::vector<Point>& points, int& batch_size, int& threads){
            std::vector<std::future<std::unordered_set<Edge>>> results;
            ThreadPool Pool(threads);
            int s = points.size();
            for(int i=0; i<s; i=i+batch_size){
                int start = i;
                int end = (i + batch_size >= s)? s: i+batch_size;
                // std::cout<<start<<" "<<end<<"\n";
                auto task = [=] {
                    std::unordered_set<Edge> visible_edges;
                    for(int j = start; j< end; j++){
                        for(auto p: visible_vertices(points[j], *G_Ptr, &points))
                        visible_edges.insert(Edge(p, points[j]));
                    }
                    return visible_edges;
                };
                results.push_back(Pool.enqueue(std::move(task)));
            }

            // progress = (status)? std::make_unique<ProgressBar>(20, results.size()): nullptr;
            for(auto &future_result: results){
                std::unordered_set<Edge> edges = future_result.get();
                for(auto edge: edges)
                VG_Ptr->add_edge(edge);
                // if (progress != nullptr) progress->update();
            }
        }

        void build(std::vector<std::vector<Point>> &polygons, std::pair<int,int> grid_resolution = {-1, -1}, int batch_size=10, int threads=1, bool status=false){
            Timer t;

            G_Ptr = std::make_unique<Graph>(polygons);
            std::vector<std::vector<Point>> empty_lvalue;
            if(grid_resolution.first <= 0) VG_Ptr = std::make_unique<Graph>(empty_lvalue);
            else VG_Ptr = std::make_unique<Graph>(empty_lvalue, grid_resolution, true);

            std::vector<Point> points(G_Ptr->points.begin(), G_Ptr->points.end());
            if(batch_size < 1) batch_size = points.size();
            // std::unique_ptr<ProgressBar> progress;

            if(threads <= 1){
                SingleThreadProcessing(points, batch_size);
            }
            else{
                MultiThreadProcessing(points, batch_size, threads);
            }
            graph_built = true;
        }

        void UpdateByADD(std::vector<std::vector<Point>> &polygons, int& batch_size, int& threads, bool tolerance=false){
            if(graph_built){
                if(VG_Ptr->record_grids){
                    std::vector<Point> new_points;
                    std::unordered_set<Point> affected_grids;

                    for(auto& polygon: polygons){
                        int s = polygon.size();
                        if (s > 1 and polygon[0] == polygon[s - 1]){
                            polygon.pop_back();
                            s -= 1;
                        }
                        int i = 0;
                        for(auto& point: polygon){
                            Point& sibling_point = polygon[(i + 1) % s];
                            Edge edge(point, sibling_point);
                            new_points.push_back(point);
                            i++;
                            BresenhamIterator bresenham(point.x, point.y, sibling_point.x, sibling_point.y);
                            while(*bresenham != bresenham.end()){
                                ++bresenham;
                                Point p((*bresenham).first, (*bresenham).second);
                                if(affected_grids.find(p) == affected_grids.end()){
                                    for(auto& visibility: VG_Ptr->section.section[p]){
                                        auto it = VG_Ptr->edges.find(visibility);
                                        if(it != VG_Ptr->edges.end()){
                                            if(tolerance || edge_intersect(visibility.p1, visibility.p2, edge)){
                                                VG_Ptr->edges.erase(it);
                                                VG_Ptr->graph[visibility.p1].erase(visibility);
                                                VG_Ptr->graph[visibility.p2].erase(visibility);
                                                VG_Ptr->section.section[p].erase(visibility);
                                            }
                                        }
                                        else{
                                            VG_Ptr->section.section[p].erase(visibility);
                                        }
                                    }
                                    affected_grids.insert(p);
                                }
                            }
                            if(s>2){
                                int pid = G_Ptr->pid +1;
                                point.polygon_id = pid;
                                sibling_point.polygon_id = pid;
                                edge.p1.polygon_id = pid;
                                edge.p2.polygon_id = pid;
                                G_Ptr->polygons[pid].insert(edge);
                            }
                            G_Ptr->add_edge(edge);
                        }
                        if(s>2) G_Ptr->pid += 1;
                    }

                    if(threads <= 1) SingleThreadProcessing(new_points, batch_size);
                    else MultiThreadProcessing(new_points, batch_size, threads); 
                }
                else{
                    std::cout<<"This operation cannot be performed since during graph building sections were not recorded\n";
                    return;
                }
            }
            else{
                std::cout<<"graph still have not been build\n";
                return;
            }
        }

        std::unordered_set<Point> find_visible_vertices(const Point& point, const Graph& graph){
            return visible_vertices(point, graph);
        }

        Graph& graph(){
            return *G_Ptr;
        }
        Graph& visgraph(){
            return *VG_Ptr;
        }
        void clear(){
            G_Ptr = nullptr;
            VG_Ptr = nullptr;
            graph_built = false;
        }


    public:
        std::unique_ptr<Graph> G_Ptr, VG_Ptr;
        bool graph_built;
};

#endif VIS_GRAPH