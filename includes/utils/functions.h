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
#include <CGAL/squared_distance_2.h>
#include <CGAL/number_utils.h> // For CGAL::to_double

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
using Face_handle = Custom_CDT::Face_handle;
using Vertex_handle = Custom_CDT::Vertex_handle;
 


//Calculate squared distance between two points
double squared_distance(const Point_2& a, const Point_2& b);

//Check if a triangle is obtuse
bool is_obtuse(const Point_2& a, const Point_2& b, const Point_2& c);

//Read JSON file
void read_json(const std::string& filename, value& jv);

//Just check if an edge is part of the additional constrains
bool is_in_constraints(const std::pair<int, int>& edge, const vector<std::pair<int, int>>& constraints);

//Just check if an edge is part of the region boundary
bool is_in_region_boundary(const std::pair<int, int>& edge, const vector<Point>& boundary);

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt, const Polygon& polygon);

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);

//Return true if approves the flip
bool can_flip(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4);

//DELETE??
void insert_points_within_boundary(CDT& cdt, const vector<Point_2>& points, const vector<Point_2>& boundary_points);
//DELETE??
bool is_point_inside_triangle(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& steiner_point);

void start_the_flips(Custom_CDT& cdt, const Polygon& polygon);

bool is_point_inside_region(const Point_2& point, const Polygon polygon);

Point find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3);

void insert_circumcenter_centroid(Custom_CDT& custom_cdt, const Polygon& polygon);

void insert_projection(Custom_CDT& custom_cdt, const Polygon polygon);

bool insert_circumcenter(Custom_CDT& custom_cdt, Custom_CDT& simulation, const Point_2& circumcenter, const Polygon& polygon);

void insert_incenter(Custom_CDT& custom_cdt, const Polygon& polygon);

void insert_bisector(Custom_CDT& custom_cdt, const Polygon& polygon);

void insert_midpoint(Custom_CDT& custom_cdt, const Polygon& polygon);

Point_2 find_incenter(const Point_2& p1, const Point_2& p2, const Point_2& p3);

bool can_insert_centroid(Custom_CDT& custom_cdt, Face_handle& triangleA, const Point_2& centroid, const Polygon& polygon);

bool can_insert_on_projection(const Point_2& projected_point, const Polygon polygon);

bool can_insert_on_midpoint(const Point_2& midpoint_point, const Polygon polygon);

Face_handle get_neighboring_face(Face_handle face, const Point_2& point1, const Point_2& point2);

Vertex_handle find_nearest_vertex(const Custom_CDT& cdt, const Point_2& query_point);

bool are_points_on_polygon(const Polygon& polygon, const Point_2& p1, const Point_2& p2);

bool circumcenter_steiner(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& circumcenter, const vector<Point_2>& points, const vector<int>& region_boundary);
bool centroid_steiner(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& centroid);
vector<Point_2> make_region_boundary(const vector<int>& region_bound_indx, const vector<Point_2>& points);


bool is_edge_in_boundary(const Point_2& p1, const Point_2& p2, const Polygon& polygon);

CGAL::Segment_2<K> find_longest_edge(const Point_2& p1, const Point_2& p2, const Point_2& p3);


