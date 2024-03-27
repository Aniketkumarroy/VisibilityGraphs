#ifndef VSISBLE_VERTICES
#define VSISBLE_VERTICES


#include<algorithm>
#include"graph.h"

const int INF = 10000;
const int CCW = 1;
const int CW = -1;
const int COLLINEAR = 0;

const int COLIN_TOLERANCE = 10;
const long T = pow(10, COLIN_TOLERANCE);
const double T2 = pow(10.0, COLIN_TOLERANCE);
const float EPS = 1e-4;

float angle(const Point& center, const Point& point);
float angle2(const Point& point_a, const Point& point_b, const Point& point_c);
float ccw(const Point& A, const Point& B, const Point& C);
bool on_segment(const Point &p, const Point &q, const Point& r);
bool edge_intersect(const Point& p1, const Point& q1, const Edge& edge);
float edge_distance(const Point &p1, const Point &p2);
Point unit_vector(const Point& c, const Point& p);
std::pair<bool, Point> intersect_point(const Point &p1, const Point &p2, const Edge& edge);

float point_edge_distance(const Point& p1, const Point& p2, const Edge& edge);
bool polygon_crossing(const Point& p1, const std::unordered_set<Edge>& poly_edges);
bool edge_in_polygon(const Point& p1, const Point& p2, const Graph& graph);
int point_in_polygon(const Point& p, const Graph& graph);
Point closest_point(const Point& p, Graph& graph, int polygon_id, double length=0.001);

class ActiveEdges{
    public:
        std::vector<Edge> active_edges;
        bool less_than(const Point& p1, const Point& p2, Edge edge1, Edge edge2);
        int index(const Point& p1, const Point& p2, const Edge& edge);

        void insert(const Point&p1, const Point&p2, Edge edge);

        void remove(const Point& p1, const Point& p2, Edge edge);

        Edge smallest();
        int size();

        Edge operator[](int index) const;
};

bool Comparator(const Point& p1, const Point& p2, const Point& point);

std::unordered_set<Point> visible_vertices(const Point& point, const Graph& graph, const std::vector<Point>* GraphPoints=nullptr, const Point* origin=nullptr, const Point* destination=nullptr, int scan=1);

#endif