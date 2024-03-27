#ifndef GRAPH
#define GRAPH

#include<iostream>
#include<string>
#include<iomanip>
#include<unordered_map>
#include<unordered_set>
#include<vector>
#include <boost/functional/hash.hpp>


/*
    This is our class wrapper for describing a point
    Point(x,y) will represent a point (x,y). this class
    properly wraps the point as a datatype of its own
    having all important functionalities.
*/
class Point{
    public:
        float x; // x co-ordinate
        float y; // y co-ordinate
        int polygon_id; // this will be required when we wil later used the collection of this points to define a polygon

        Point(float x, float y, int polygon_id=-1) : x(x), y(y), polygon_id(polygon_id) {}
        bool operator==(const Point& point) const {
            return x == point.x && y == point.y;
        }
        bool operator!=(const Point& point) const {
            return !(*this == point);
        }
        bool operator<(const Point& point) const {
            std::hash<float> hasher;
            std::size_t hash1 = hasher(this->x) ^ hasher(this->y);
            std::size_t hash2 = hasher(point.x) ^ hasher(point.y);
            return hash1 < hash2;
        }
        std::string toString() const {
            std::ostringstream oss;
            oss << "Point(" << std::fixed << std::setprecision(2) << x << ", " << y << ", " << polygon_id << ")";
            return oss.str();
        }
};

/*
    This is our class wrapper for describing a Edge.
    an edge is described by two points. this class
    properly wraps the edge as a datatype of its own
    having all important functionalities.
*/
class Edge{
    public:
        Point p1;
        Point p2;
        Edge(Point point1, Point point2) : p1(point1), p2(point2) {}
        Point get_adjacent(const Point& point) const {
            // given a point of the edge, it returns the other point, or its sibling point
            if (point == p1) return p2;
            return p1;
        }
        bool contains(const Point& point) const {
            // checks whether the edge contains the given point
            return (point == p1) || (point == p2);
        }
        bool operator==(const Edge& edge) const {
            if (this->p1 == edge.p1 and this->p2 == edge.p2)
            return true;
            else if (this->p2 == edge.p1 and this->p1 == edge.p2)
            return true;
            else return false;
        }
        bool operator!=(const Edge& edge) const {
            return !(*this==edge);
        }
        std::string toString() const {
            std::ostringstream oss;
            oss << "Edge(" << p1.toString()<<", "<<p2.toString()<<")";
            return oss.str();
        }
};

namespace std {
    template<>
    struct hash<Point> {
        std::size_t operator()(const Point& point) const {
            // Define your custom hashing logic here
            /*
                note that we can start with a XOR hashing function like
                return std::hash<float>()(point.x) ^ std::hash<float>()(point.y);
                but this will produce collisions for points such as Point(a, b) and Point(b, a)
                since Point(a,b) != Point(b,a), we want hash(Point(a, b)) != hash(Point(b,a))
            */

            // https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
            std::size_t seed = 0;
            boost::hash_combine(seed, boost::hash_value(point.x));
            boost::hash_combine(seed, boost::hash_value(point.y));
            return seed;
        }
    };
    template<>
    struct hash<Edge> {
        std::size_t operator()(const Edge& edge) const {
            std::hash<Point> hasher;
            // since a Edge(p1, p2) = Edge(p2, p1) we will use XOR for hashing
            return hasher(edge.p1) ^ hasher(edge.p2);
        }
    };
};

class BresenhamIterator{
    private:
        int x_start, x_end, y_start, y_end;
        int _dx, _dy, _xx, _yy, _xy, _yx, _x, _y, _xsign, _ysign, _D;
        std::pair<int, int> P;
        std::pair<int, int> Begin;
        std::pair<int, int> End;
    public:
        BresenhamIterator(int x_start, int y_start, int x_end, int y_end){
            this->x_start = x_start;
            this->y_start = y_start;
            this->x_end =x_end;
            this->y_end = y_end;
            Begin = {x_start, y_start};
            End = {x_end, y_end};

            _dx = x_end - x_start;
            _dy = y_end - y_start;
            _xsign = (_dx > 0)?1:-1;
            _ysign = (_dy > 0)?1:-1;
            _dx = abs(_dx);
            _dy = abs(_dy);
            if (_dx > _dy){
                _xx = _xsign;
                _xy = 0;
                _yx = 0;
                _yy = _ysign;
            }else{
                std::swap(_dx, _dy);
                _xx = 0;
                _xy = _ysign;
                _yx = _xsign;
                _yy = 0;
            }
            _D = 2*_dy - _dx;
            _y = 0;
            _x = 0;
            P.first = x_start;
            P.second = y_start;
        }
    
    std::pair<int, int> operator*(){
        return P;
    }

    void operator++(){
        if (_x <= _dx){
            P.first = (x_start + _x*_xx + _y*_yx);
            P.second = (y_start + _x*_xy + _y*_yy);
            if (_D >= 0){
                _y += 1;
                _D -= 2*_dx;
                }
            _D += 2*_dy;
            _x += 1;
        }
        else
        P = End;
    }

    std::pair<int, int> begin(){
        return Begin;
    }
    std::pair<int, int> end(){
        return End;
    }
};

/*
    This is a datastructure for storing the grid as a key and 
    the set of edges that passess from it as a value. a grid is defined by
    its one corner point and its height and width
*/

class Section{
    public:
        std::unordered_map<Point, std::unordered_set<Edge>> section; // internal datastructure for storing
        std::pair<int, int> resolution; // grid resolution: (grid_width, grid_height)
        void set_grid_resolution(std::pair<int, int> grid_resolution){
            resolution = grid_resolution;
        }
        std::unordered_set<Edge> operator[](const Point& point){
            // overloading the index operator to get the edges for the given point
            // if the point is not found, a empty set will be returned
            return section[point];
        }
        void add(const Edge& edge){
            int x0 = (int)(edge.p1.x/resolution.first);
            int y0 = (int)(edge.p1.y/resolution.second);
            int x1 = (int)(edge.p2.x/resolution.first);
            int y1 = (int)(edge.p2.y/resolution.second);

            BresenhamIterator bresenham(x0, y0, x1, y1);
            while(*bresenham != bresenham.end()){
                ++bresenham;
                section[Point((*bresenham).first, (*bresenham).second)].insert(edge);
            }
        }
};

class Graph{
    public:
        std::unordered_map<Point, std::unordered_set<Edge>> graph;
        std::unordered_set<Edge> edges;
        std::unordered_set<Point> points;
        std::unordered_map<int, std::unordered_set<Edge>> polygons;
        Section section;
        std::pair<int, int> grid = {-1, -1};
        bool record_grids = false;
        int pid = 0;
        Graph(std::vector<std::vector<Point>> &polygons, std::pair<int, int> grid_resolution={1, 1}, bool record_section=false){
            grid = grid_resolution;
            record_grids = record_section;
            if (record_grids)
            section.set_grid_resolution(grid);

            for(auto &polygon: polygons){
                int s = polygon.size();
                if (s > 1 and polygon[0] == polygon[s - 1]){
                    polygon.pop_back();
                    s -= 1;
                }
                int i = 0;
                for(auto &point: polygon){
                    Point& sibling_point = polygon[(i + 1) % s];
                    Edge edge(point, sibling_point);
                    if (s > 2){
                        point.polygon_id = pid;
                        sibling_point.polygon_id = pid;
                        edge.p1.polygon_id = pid;
                        edge.p2.polygon_id = pid;
                        this->polygons[pid].insert(edge);
                    }
                    this->add_edge(edge);
                    i++;
                }
                if (s > 2)
                pid += 1;
            }
        }
        std::vector<Point> get_adjacent_points(const Point& point) const {
            std::vector<Point> v;
            try{
                for(const auto &edge: graph.at(point))
                v.push_back(edge.get_adjacent(point));
            }catch (const std::out_of_range& e) {
                // Key not found, return an empty set
                return v;
            }
            return v;
        }
        void add_edge(const Edge &edge){
            graph[edge.p1].insert(edge);
            graph[edge.p2].insert(edge);
            points.insert(edge.p1);
            points.insert(edge.p2);
            edges.insert(edge);
            if (record_grids)
                section.add(edge);
        }
        bool contains(const Point &point) const{
            if (points.find(point) != points.end())
            return true;
            return false;
        }
        bool contains(const Edge &edge) const{
            if(edges.find(edge) != edges.end())
            return true;
            return false;
        }
        std::unordered_set<Edge> operator[](const Point& point) const {
            try{
                // Attempt to access the set directly using at()
                return graph.at(point);
            }catch (const std::out_of_range& e) {
                // Key not found, return an empty set
                return std::unordered_set<Edge>{};
            }
        }
        std::unordered_set<Edge> operator[](int polygon_id) const {
            try{
                // Attempt to access the set directly using at()
                return polygons.at(polygon_id);
            }catch (const std::out_of_range& e) {
                // Key not found, return an empty set
                return std::unordered_set<Edge>{};
            }
        }

        std::string toString(){
            std::string exp = "";
            for(const auto &point: points){
                exp += point.toString() + ":\n";
                for(const auto &edge: graph[point]){
                    exp += "\t" + edge.toString() + "\n";
                }
            }
            return exp;
        }
};

#endif GRAPH