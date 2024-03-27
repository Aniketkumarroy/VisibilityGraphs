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

        Point(float x, float y, int polygon_id=-1);
        bool operator==(const Point& point) const;
        bool operator!=(const Point& point) const;
        bool operator<(const Point& point) const;
        std::string toString() const;
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
        Edge(Point point1, Point point2);
        Point get_adjacent(const Point& point) const;
        bool contains(const Point& point) const;
        bool operator==(const Edge& edge) const;
        bool operator!=(const Edge& edge) const;
        std::string toString() const;
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
        BresenhamIterator(int x_start, int y_start, int x_end, int y_end);
    
        std::pair<int, int> operator*();

        void operator++();

        std::pair<int, int> begin();
        std::pair<int, int> end();
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
        void set_grid_resolution(std::pair<int, int> grid_resolution);
        std::unordered_set<Edge> operator[](const Point& point);
        void add(const Edge& edge);
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
        Graph(std::vector<std::vector<Point>> &polygons, std::pair<int, int> grid_resolution={1, 1}, bool record_section=false);
        std::vector<Point> get_adjacent_points(const Point& point) const;
        void add_edge(const Edge &edge);
        bool contains(const Point &point) const;
        bool contains(const Edge &edge) const;
        std::unordered_set<Edge> operator[](const Point& point) const;
        std::unordered_set<Edge> operator[](int polygon_id) const;

        std::string toString();
};

#endif