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

float angle(const Point& center, const Point& point){
    float dx = point.x - center.x;
    float dy = point.y - center.y;
    if (dx == 0){
        if (dy<0) return M_PI*3/2;
        else return M_PI/2;
    }else if (dy==0){
        if (dx<0) return M_PI;
        else return 0;
    }
    if (dx<0) return M_PI + atan(dy/dx);
    if (dy<0) return 2*M_PI + atan(dy/dx);
    return atan(dy/dx);
}
float angle2(const Point& point_a, const Point& point_b, const Point& point_c){
    float a = pow((point_c.x - point_b.x), 2) + pow((point_c.y - point_b.y), 2);
    float b = pow((point_c.x - point_a.x), 2) + pow((point_c.y - point_a.y), 2);
    float c = pow((point_b.x - point_a.x), 2) + pow((point_b.y - point_a.y), 2);
    float cos_value = (a + c - b) / (2 * sqrt(a) * sqrt(c));
    return acos((cos_value*T)/T2);
}
float ccw(const Point& A, const Point& B, const Point& C){
    float area = (((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x))*T)/T2;
    if (area>0) return CCW;
    if (area<0) return CW;
    return COLLINEAR;
}
bool on_segment(const Point &p, const Point &q, const Point& r){
    // Given three colinear points p, q, r, the function checks if point q lies on line segment 'pr'.
    if (q.x <= std::max(p.x, r.x) && q.x >= std::min(p.x, r.x)) {
        if (q.y <= std::max(p.y, r.y) && q.y >= std::min(p.y, r.y))
        return true;
    }
    return false;
}
bool edge_intersect(const Point& p1, const Point& q1, const Edge& edge){
    // http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    // check whether line segment from p1 and p2 intersect edge or not
    Point p2 = edge.p1;
    Point q2 = edge.p2;
    int o1 = ccw(p1, q1, p2);
    int o2 = ccw(p1, q1, q2);
    int o3 = ccw(p2, q2, p1);
    int o4 = ccw(p2, q2, q1);

    if (o1 != o2 && o3 != o4)
    return true;
    // p1, q1 and p2 are colinear and p2 lies on segment p1q1
    if (o1 == COLLINEAR && on_segment(p1, p2, q1))
    return true;
    // p1, q1 and p2 are colinear and q2 lies on segment p1q1
    if (o2 == COLLINEAR && on_segment(p1, q2, q1))
    return true;
    // p2, q2 and p1 are colinear and p1 lies on segment p2q2
    if (o3 == COLLINEAR && on_segment(p2, p1, q2))
    return true;
    // p2, q2 and q1 are colinear and q1 lies on segment p2q2
    if (o4 == COLLINEAR && on_segment(p2, q1, q2))
    return true;
    return false;
}
float edge_distance(const Point &p1, const Point &p2){
    return sqrt(pow((p1.x-p2.x), 2) + pow((p1.y-p2.y), 2));
}
Point unit_vector(const Point& c, const Point& p){
    float magnitute = edge_distance(c, p);
    return Point((p.x-c.x)/magnitute, (p.y-c.y)/magnitute);
}
std::pair<bool, Point> intersect_point(const Point &p1, const Point &p2, const Edge& edge){
    if (edge.contains(p1)) return std::pair<bool, Point> {true, p1};
    if (edge.contains(p2)) return std::pair<bool, Point> {true, p2};
    if (edge.p1.x == edge.p2.x){
        if (p1.x == p2.x)
        return std::pair<bool, Point> {false, Point(-1, -1)};
        float pslope = (p1.y - p2.y) / (p1.x - p2.x);
        float intersect_x = edge.p1.x;
        float intersect_y = pslope * (intersect_x - p1.x) + p1.y;
        return std::pair<bool, Point> {true, Point(intersect_x, intersect_y)};
    }
    if (p1.x == p2.x){
        float eslope = (edge.p1.y - edge.p2.y) / (edge.p1.x - edge.p2.x);
        float intersect_x = p1.x;
        float intersect_y = eslope * (intersect_x - edge.p1.x) + edge.p1.y;
        return std::pair<bool, Point> {true, Point(intersect_x, intersect_y)};
    }
    float pslope = (p1.y - p2.y) / (p1.x - p2.x);
    float eslope = (edge.p1.y - edge.p2.y) / (edge.p1.x - edge.p2.x);
    if (eslope == pslope)
    return std::pair<bool, Point> {false, Point(-1, -1)};
    float intersect_x = (eslope * edge.p1.x - pslope * p1.x + p1.y - edge.p1.y) / (eslope - pslope);
    float intersect_y = eslope * (intersect_x - edge.p1.x) + edge.p1.y;
    return std::pair<bool, Point> {true, Point(intersect_x, intersect_y)};
}

float point_edge_distance(const Point& p1, const Point& p2, const Edge& edge){
    std::pair<bool, Point> ip = intersect_point(p1, p2, edge);
    if(!ip.first)
    return 0;
    return edge_distance(p1, ip.second);
}
bool polygon_crossing(const Point& p1, const std::unordered_set<Edge>& poly_edges){
    // Returns True if Point p1 is internal to the polygon. The polygon is
    // defined by the Edges in poly_edges. Uses crossings algorithm and takes into
    // account edges that are collinear to p1.
    // https://www.geeksforgeeks.org/how-to-check-if-a-given-point-lies-inside-a-polygon/
    Point p2(INF, p1.y);
    int intersect_count = 0;
    for(const auto &edge: poly_edges){
        if (p1.y < edge.p1.y && p1.y < edge.p2.y) continue; // if the points lies entirely above the edge
        if (p1.y > edge.p1.y && p1.y > edge.p2.y) continue; // if the points lies entirely above the edge
        if (p1.x > edge.p1.x && p1.x > edge.p2.x) continue; // if the points lies entirely right the edge

        bool edge_p1_collinear = (ccw(p1, edge.p1, p2) == COLLINEAR);
        bool edge_p2_collinear = (ccw(p1, edge.p2, p2) == COLLINEAR);
        if (edge_p1_collinear && edge_p2_collinear) continue; // if the edge is coincidenet to p1->p2
        if (edge_p1_collinear || edge_p2_collinear){
            Point collinear_point = (edge_p1_collinear) ? edge.p1 : edge.p2;
            if (edge.get_adjacent(collinear_point).y > p1.y)
            intersect_count += 1;
        }else if (edge_intersect(p1, p2, edge))
            intersect_count += 1;
    }
    if (intersect_count % 2 == 0)
    return false;
    return true;
}
bool edge_in_polygon(const Point& p1, const Point& p2, const Graph& graph){
    //Return true if the edge from p1 to p2 is completely interior to any polygon in graph.
    if (p1.polygon_id != p2.polygon_id) return false;
    if (p1.polygon_id == -1 or p2.polygon_id == -1) return false;
    Point mid_point((p1.x + p2.x) / 2, (p1.y + p2.y) / 2);
    return polygon_crossing(mid_point, graph[p1.polygon_id]);
}
int point_in_polygon(const Point& p, const Graph& graph){
    // Return true if the point p is interior to any polygon in graph.
    for(const auto &polygon: graph.polygons){
        if (polygon_crossing(p, polygon.second))
        return polygon.first;
    }
    return -1;
}
Point closest_point(const Point& p, Graph& graph, int polygon_id, double length=0.001){
    // Assumes p is interior to the polygon with polygon_id. Returns the
    // closest point c outside the polygon to p, where the distance from c to
    // the intersect point from p to the edge of the polygon is length.

    Point close_point(0.0, 0.0);
    Edge close_edge(close_point, close_point);
    float close_dist;
    // Finds point closest to p, but on a edge of the polygon.
    // Solution from http://stackoverflow.com/a/6177788/4896361
    int i = 0;
    for(const auto &e: graph[polygon_id]){
        float num = ((p.x-e.p1.x)*(e.p2.x-e.p1.x) + (p.y-e.p1.y)*(e.p2.y-e.p1.y));
        float denom = (pow((e.p2.x - e.p1.x),2) + pow((e.p2.y - e.p1.y),2));
        float u = num/denom;
        Point pu(e.p1.x + u*(e.p2.x - e.p1.x), e.p1.y + u*(e.p2.y- e.p1.y));
        Point pc = pu;
        if (u < 0) pc = e.p1;
        else if (u > 1) pc = e.p2;
        float d = edge_distance(p, pc);
        if (i == 0 || d < close_dist){
            close_dist = d;
            close_point = pc;
            close_edge = e;
        }
        i++;
    }
    if (close_edge.contains(close_point)){
        Point c = (close_point == close_edge.p1) ? close_edge.p1 : close_edge.p2;
        std::vector<Point> adjacent_points = graph.get_adjacent_points(c);
        Edge e1(c, adjacent_points[0]);
        Edge e2(c, adjacent_points[1]);
        Point v1 = unit_vector(c, e1.get_adjacent(c));
        Point v2 = unit_vector(c, e2.get_adjacent(c));
        Point vsum = unit_vector(Point(0, 0), Point(v1.x + v2.x, v1.y + v2.y));
        Point close1(c.x + (vsum.x * length), c.y + (vsum.y * length));
        Point close2(c.x - (vsum.x * length), c.y - (vsum.y * length));
        if (point_in_polygon(close1, graph) == -1)
            return close1;
        return close2;
    }else{
        Point v = unit_vector(p, close_point);
        return Point(close_point.x + v.x*length, close_point.y + v.y*length);
    }
}

class ActiveEdges{
    public:
        std::vector<Edge> active_edges;
        bool less_than(const Point& p1, const Point& p2, Edge edge1, Edge edge2){
            if (edge1 == edge2)
                return false;
            if (!edge_intersect(p1, p2, edge2))
                return true;
            float edge1_dist = point_edge_distance(p1, p2, edge1);
            float edge2_dist = point_edge_distance(p1, p2, edge2);
            if (edge1_dist > edge2_dist)
                return false;
            if (edge1_dist < edge2_dist)
                return true;
            if (abs(edge1_dist - edge2_dist) <= 0.01){ // comparing float by having a tolerance of 0.01
                Point same_point = edge1.p1;
                if (edge2.contains(edge1.p1))
                    same_point = edge1.p1;
                else
                    same_point = edge1.p2;
                float angle_edge1 = angle2(p1, p2, edge1.get_adjacent(same_point));
                float angle_edge2 = angle2(p1, p2, edge2.get_adjacent(same_point));
                if (angle_edge1 < angle_edge2)
                    return true;
                return false;
            }
            return false;
        }
        int index(const Point& p1, const Point& p2, const Edge& edge){
            int lo = 0;
            int hi = active_edges.size();
            while (lo < hi){
                int mid = (lo+hi)/2;
                if (this->less_than(p1, p2, edge, active_edges[mid]))
                    hi = mid;
                else
                    lo = mid + 1;
            }
            return lo;
        }

        void insert(const Point&p1, const Point&p2, Edge edge){
            int iDx = this->index(p1, p2, edge);
            active_edges.insert(active_edges.begin() + iDx, edge);
        }

        void remove(const Point& p1, const Point& p2, Edge edge){
            int iDx = this->index(p1, p2, edge) - 1;
            if (active_edges[iDx] == edge)
                active_edges.erase(active_edges.begin() + iDx);
        }

        Edge smallest(){
            return active_edges[0];
        }

        int size(){
            return active_edges.size();
        }

        Edge operator[](int index) const {
            return active_edges[index];
        }
};

bool Comparator(const Point& p1, const Point& p2, const Point& point) {
    double angle1 = angle(point, p1);
    double angle2 = angle(point, p2);

    if (angle1 != angle2) {
        return angle1 < angle2;
    }

    double distance1 = edge_distance(point, p1);
    double distance2 = edge_distance(point, p2);

    return distance1 < distance2;
}

std::unordered_set<Point> visible_vertices(const Point& point, const Graph& graph, const std::vector<Point>* GraphPoints=nullptr, const Point* origin=nullptr, const Point* destination=nullptr, int scan=1){
    const std::unordered_set<Edge>& edges = graph.edges;
    std::vector<Point> points = (GraphPoints != nullptr) ? *GraphPoints : std::vector<Point>(graph.points.begin(), graph.points.end());
    if (origin != nullptr) points.push_back(*origin);
    if (destination != nullptr) points.push_back(*destination);
    std::sort(points.begin(), points.end(), [point](Point p1, Point p2){return Comparator(p1, p2, point);});

    // Initialize open_edges with any intersecting edges on the half line from
    // point along the positive x-axis
    ActiveEdges active;
    Point point_inf(INF, point.y);
    for(const auto &edge: edges){
        if (edge.contains(point)) continue;
        if (edge_intersect(point, point_inf, edge)){
            if (on_segment(point, edge.p1, point_inf)) continue;
            if (on_segment(point, edge.p2, point_inf)) continue;
            active.insert(point, point_inf, edge);
        }
    }

    std::unordered_set<Point> visible;
    Point prev(-1, -1, -2);
    bool prev_visible = false;
    bool prev_available = false;
    for(const auto &p: points){
        if (p == point) continue;
        if (scan == 0 and angle(point, p) > M_PI) break;

        // Update open_edges - remove clock wise edges incident on p
        if (active.size() > 0){
            for(const auto &edge: graph[p]){
                if (ccw(point, p, edge.get_adjacent(p)) == CW)
                    active.remove(point, p, edge);
            }
        }

        // Check if p is visible from point
        bool is_visible = false;
        // ...Non-collinear points
        if (!prev_available || ccw(point, prev, p) != COLLINEAR || !on_segment(point, prev, p)){
            if (active.size() == 0)
                is_visible = true;
            else if (!edge_intersect(point, p, active.smallest()))
                is_visible = true;
        }
        // ...For collinear points, if previous point was not visible, p is not
        else if (!prev_visible)
            is_visible = false;
        // ...For collinear points, if previous point was visible, need to check
        // that the edge from prev to p does not intersect any open edge.
        else{
            is_visible = true;
            for (const auto &edge: active.active_edges){
                if (!edge.contains(prev) && edge_intersect(prev, p, edge)){
                    is_visible = false;
                    break;
                }
            }
            if (is_visible && edge_in_polygon(prev, p, graph))
                    is_visible = false;
        }
        // Check if the visible edge is interior to its polygon
        if (is_visible){
            std::vector<Point> adjacent_points = graph.get_adjacent_points(point);
            if(find(adjacent_points.begin(), adjacent_points.end(), p) == adjacent_points.end())
                is_visible = !edge_in_polygon(point, p, graph);
        }

        if (is_visible) visible.insert(p);

        // Update open_edges - Add counter clock wise edges incident on p
        for(const auto &edge: graph[p]){
            if (!edge.contains(point) && ccw(point, p, edge.get_adjacent(p)) == CCW)
                active.insert(point, p, edge);
        }

        prev = p;
        prev_available = true;
        prev_visible = is_visible;
    }
    return visible;
}

#endif VSISBLE_VERTICES