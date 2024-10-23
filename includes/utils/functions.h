#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <iostream>
#include <vector>
#include <boost/json.hpp>
#include <fstream>
#include <string>
#include <cmath>
#include <CGAL/Polygon_2.h>
//#include <CGAL/bounded_side_2.h> // For bounded_side_2
//#include <CGAL/circumcenter.h>
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"

using namespace boost::json;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Polygon = CGAL::Polygon_2<K>;
//typedef CGAL::Polygon_2<K> Polygon;


//Calculate squared distance between two points
double squared_distance(const Point& a, const Point& b);

//Check if a triangle is obtuse
bool is_obtuse(const Point& a, const Point& b, const Point& c);

//Read JSON file
void read_json(const std::string& filename, value& jv);

//Just check if an edge is part of the additional constrains
bool is_in_constraints(const std::pair<int, int>& edge, const std::vector<std::pair<int, int>>& constraints);

//Just check if an edge is part of the region boundary
bool is_in_region_boundary(const std::pair<int, int>& edge, const std::vector<Point>& boundary);

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt);

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point& p1, const Point& p2, const Point& p3, const Point& p4);

//Return true if approves the flip
bool can_flip(const Point& p1, const Point& p2, const Point& p3, const Point& p4);

void insert_points_within_boundary(CDT& cdt, const std::vector<Point>& points, const std::vector<Point>& boundary_points);

//Ιnsert Steiner points at circumcenters
void insert_steiner_points(Custom_CDT& custom_cdt, const std::vector<Point>& points, const std::vector<int>& region_boundary);

bool is_point_inside_triangle(const Point& p1, const Point& p2, const Point& p3, const Point& steiner_point);

void start_the_flips(Custom_CDT& cdt,const std::vector<Point>& points,const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary);

bool is_point_inside_region(const Point& point, const std::vector<int>& region_boundary_indx, const std::vector<Point>& points);

bool circumcenter_steiner(const Point& p1, const Point& p2, const Point& p3, const Point& circumcenter, const std::vector<Point>& points, const std::vector<int>& region_boundary);
bool centroid_steiner(const Point& p1, const Point& p2, const Point& p3, const Point& centroid);
std::vector<Point> make_region_boundary(const std::vector<int>& region_bound_indx, const std::vector<Point>& points);