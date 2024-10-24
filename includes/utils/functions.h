#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <boost/json.hpp>
#include <fstream>
#include <string>
#include <cmath>
#include <CGAL/Polygon_2.h>
#include <CGAL/number_utils.h>
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Line_2.h>

using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
//using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Polygon = CGAL::Polygon_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;
 


//Calculate squared distance between two points
double squared_distance(const Point_2& a, const Point_2& b);

//Check if a triangle is obtuse
bool is_obtuse(const Point_2& a, const Point_2& b, const Point_2& c);

//Read JSON file
void read_json(const std::string& filename, value& jv);

//Just check if an edge is part of the additional constrains
bool is_in_constraints(const std::pair<int, int>& edge, const std::vector<std::pair<int, int>>& constraints);

//Just check if an edge is part of the region boundary
bool is_in_region_boundary(const std::pair<int, int>& edge, const std::vector<Point>& boundary);

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt);

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);

//Return true if approves the flip
bool can_flip(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);

void insert_points_within_boundary(CDT& cdt, const std::vector<Point_2>& points, const std::vector<Point_2>& boundary_points);

bool is_point_inside_triangle(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& steiner_point, const Polygon polygon);

void start_the_flips(Custom_CDT& cdt,const std::vector<Point_2>& points,const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary);

bool is_point_inside_region(const Point_2& point, const Polygon polygon);

Point find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3);

void insert_steiner_points(Custom_CDT& custom_cdt, const std::vector<Point_2>& points, const std::vector<int>& region_boundary, const Polygon polygon);

bool can_insert_circumcenter(Custom_CDT& custom_cdt, Custom_CDT::Face_handle& triangleA, const Point_2& circumcenter);

bool can_insert_centroid(Custom_CDT& custom_cdt, Custom_CDT::Face_handle& triangleA, const Point_2& centroid);

bool can_insert_on_projection(const Point_2& projected_point, const Polygon polygon);

Custom_CDT::Face_handle get_neighboring_face(Custom_CDT::Face_handle face, const Point_2& point1, const Point_2& point2);

bool circumcenter_steiner(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& circumcenter, const std::vector<Point_2>& points, const std::vector<int>& region_boundary);
bool centroid_steiner(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& centroid);
std::vector<Point_2> make_region_boundary(const std::vector<int>& region_bound_indx, const std::vector<Point_2>& points);