#include "functions.h"
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"
#include <CGAL/Polygon_2_algorithms.h>

using namespace boost::json;
using namespace std;
//using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
typedef CGAL::Polygon_2<K> Polygon;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;
using Face_handle = Custom_CDT::Face_handle;
using Vertex_handle = Custom_CDT::Vertex_handle;

//Calculate squared distance between two points
double squared_distance(const Point_2& a, const Point_2& b) {
    return CGAL::to_double(CGAL::squared_distance(a, b));
}

//Check if a triangle is obtuse
bool is_obtuse(const Point_2& a, const Point_2& b, const Point_2& c) {
    double ab2 = CGAL::to_double(squared_distance(a, b));
    double ac2 = CGAL::to_double(squared_distance(a, c));
    double bc2 = CGAL::to_double(squared_distance(b, c));
    return (ab2 + ac2 < bc2) || (ab2 + bc2 < ac2) || (ac2 + bc2 < ab2);
}

//Read JSON file
void read_json(const std::string& filename, value& jv) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error opening file: " << filename <<endl;
        return;
    }
    std::string json_str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    jv = parse(json_str);
}

//Just check if an edge is part of the additional constrains
bool is_in_constraints(const std::pair<int, int>& edge, const vector<std::pair<int, int>>& constraints) {
    return std::find(constraints.begin(), constraints.end(), edge) != constraints.end() ||
           std::find(constraints.begin(), constraints.end(), std::make_pair(edge.second, edge.first)) != constraints.end();
           //checks for the edge in reverse order as well

}

//Just check if an edge is part of the region boundary
bool is_in_region_boundary(const std::pair<int, int>& edge, const vector<int>& boundary) {
    for (int i = 0; i < boundary.size(); ++i) {
        int idx1 = boundary[i];
        int idx2 = boundary[(i + 1) % boundary.size()];
        if ((edge.first == idx1 && edge.second == idx2) || (edge.first == idx2 && edge.second == idx1)) {
            return true;
        }
    }
    return false;
}

//Just count the number of obtuses triangles in a cdt
int count_obtuse_triangles(CDT& cdt, const Polygon& polygon) {
    int obtuse_count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();

        if (is_obtuse(p1, p2, p3) && 
            polygon.bounded_side(CGAL::midpoint(p1, p2)) != CGAL::ON_UNBOUNDED_SIDE && 
            polygon.bounded_side(CGAL::midpoint(p1, p3)) != CGAL::ON_UNBOUNDED_SIDE && 
            polygon.bounded_side(CGAL::midpoint(p2, p3)) != CGAL::ON_UNBOUNDED_SIDE) {
            obtuse_count++;
        }
    }
    return obtuse_count;
}

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4) {
    vector<Point_2> points = {p1, p2, p3, p4};

    //Check if the polygon is convex
    return CGAL::is_convex_2(points.begin(), points.end(), K()); 
}

//Return true if approves the flip
bool can_flip(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4){
    //Check if the quadrilateral (p2, p4, p3, p1) is convex
    if (is_convex(p1, p2, p3, p4)) {
            // Count obtuse angles before the flip
            int obtuse_before_cnt = 0;
            if(is_obtuse(p1, p2, p3)) obtuse_before_cnt++;
            if(is_obtuse(p1, p3, p4)) obtuse_before_cnt++;

            //Count obtuse angles after the potential flip
            int obtuse_after_cnt = 0;
            if(is_obtuse(p1, p2, p4)) obtuse_after_cnt++;
            if(is_obtuse(p2, p3, p4)) obtuse_after_cnt++;

            //Check collinearity of all possible triplets among p1, p2, p3, p4
            if (CGAL::orientation(p1, p2, p3) == CGAL::COLLINEAR ||
                CGAL::orientation(p1, p2, p4) == CGAL::COLLINEAR ||
                CGAL::orientation(p1, p3, p4) == CGAL::COLLINEAR ||
                CGAL::orientation(p2, p3, p4) == CGAL::COLLINEAR) {
                //Skip the flip if any three points are collinear
                return false;
            }

            //If its worth flipping, do it!
            if(obtuse_after_cnt < obtuse_before_cnt) return true;
            else return false;
    }
    //Case quadrilateral is not convex
    return false;      
}

//Midpoint Insertion:
//Finds the longest edge of the obtuse triangle and calculates its midpoint.
//Simulates inserting the midpoint to see if it reduces obtuse angles. If successful, inserts it into the actual CDT.
void insert_midpoint(Custom_CDT& custom_cdt, const Polygon& polygon) {
    bool progress = true;
    while (progress) {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            if (is_obtuse(p1, p2, p3)) {
                Custom_CDT simulation = custom_cdt;
                CGAL::Segment_2 longest_edge = find_longest_edge(p1, p2, p3);
                Point_2 midpoint = CGAL::midpoint(longest_edge.source(), longest_edge.target());

                if (can_insert_on_midpoint(midpoint, polygon)) {
                    int obtuses_before = count_obtuse_triangles(simulation, polygon);
                    simulation.insert_no_flip(midpoint);
                    start_the_flips(simulation, polygon);
                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)) {
                        custom_cdt.insert_no_flip(midpoint);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        cout << "Midpoint inserted to reduce obtuse angles.\n";
                        break;
                    }
                }
            }
        }
    }
}

bool can_insert_on_midpoint(const Point_2& midpoint_point, const Polygon polygon) {
    //Check if you projected point is on side (boundary)
    if (is_point_inside_region(midpoint_point, polygon)) return true;
    else return false;
}


void insert_incenter(Custom_CDT& custom_cdt, const Polygon& polygon) {
    bool progress = true;
    while (progress) {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            if (is_obtuse(p1, p2, p3)) {
                Custom_CDT simulation = custom_cdt;
                Point_2 incenter = find_incenter(p1, p2, p3);

                if (can_insert_on_midpoint(incenter, polygon)) {
                    int obtuses_before = count_obtuse_triangles(simulation, polygon);
                    simulation.insert_no_flip(incenter);
                    start_the_flips(simulation, polygon);
                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)) {
                        custom_cdt.insert_no_flip(incenter);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        std::cout << "Incenter inserted to reduce obtuse angles.\n";
                        break;
                    }
                }
            }
        }
    }
}


Point_2 find_incenter(const Point_2& p1, const Point_2& p2, const Point_2& p3) {
    // Calculate the lengths of the sides
    double a = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3)));
    double b = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p3)));
    double c = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));

    // Calculate the coordinates of the incenter using CGAL types
    auto incenter_x = (a * CGAL::to_double(p1.x()) + b * CGAL::to_double(p2.x()) + c * CGAL::to_double(p3.x())) / (a + b + c);
    auto incenter_y = (a * CGAL::to_double(p1.y()) + b * CGAL::to_double(p2.y()) + c * CGAL::to_double(p3.y())) / (a + b + c);

    // Return the incenter as a Point_2 object
    return Point_2(incenter_x, incenter_y);
}

void insert_bisector(Custom_CDT& custom_cdt, const Polygon& polygon) {
    bool progress = true;
    while (progress) {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            if (is_obtuse(p1, p2, p3)) {
                Custom_CDT simulation = custom_cdt;
                Point_2 obtuse_vertex = find_obtuse_vertex(p1, p2, p3);
                Point_2 opposite1, opposite2;

                // Identify the opposite vertices based on the obtuse vertex
                if (obtuse_vertex == p1) {
                    opposite1 = p2;
                    opposite2 = p3;
                } else if (obtuse_vertex == p2) {
                    opposite1 = p1;
                    opposite2 = p3;
                } else {
                    opposite1 = p1;
                    opposite2 = p2;
                }

                // Calculate the centroid of the triangle formed by the three points
                Point_2 centroid_opposites = CGAL::centroid(obtuse_vertex, opposite1, opposite2);
                // Calculate the midpoint between the obtuse vertex and the centroid
                Point_2 bisector_midpoint = CGAL::midpoint(obtuse_vertex, centroid_opposites);

                if (can_insert_on_projection(bisector_midpoint, polygon)) {
                    int obtuses_before = count_obtuse_triangles(simulation, polygon);
                    simulation.insert_no_flip(bisector_midpoint);
                    start_the_flips(simulation, polygon);
                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)) {
                        custom_cdt.insert_no_flip(bisector_midpoint);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        std::cout << "Bisector midpoint inserted to reduce obtuse angles.\n";
                        break;
                    }
                }
            }
        }
    }
}




static int num = 0;
static int obt = 0;
void insert_projection(Custom_CDT& custom_cdt, const Polygon polygon){
    //Projection case
    bool progress = true;
    while(progress){
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face){
            num++;
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();
            Point_2 opposite1, opposite2;

            //Copy the cdt for simulation
            Custom_CDT simulation = custom_cdt;
            if(is_obtuse(p1, p2, p3)){
                obt++;
                Point_2 obtuse_angle_vertex = find_obtuse_vertex(p1, p2, p3);            
                //Face_handle triangleA = face;
                //Find the obtuse point, and the 2 opposites
                if (obtuse_angle_vertex == p1) {
                    opposite1 = p2;
                    opposite2 = p3;
                } else if (obtuse_angle_vertex == p2) {
                    opposite1 = p1;
                    opposite2 = p3;
                } else {
                    opposite1 = p1;
                    opposite2 = p2;
                }
                Line_2 line(opposite1, opposite2);
                Point_2 projected_point = line.projection(obtuse_angle_vertex);
                bool insert_projection = can_insert_on_projection(projected_point, polygon);
                if(insert_projection){
                    /*Simulate insertion of Projection*/
                    int obtuses_before = count_obtuse_triangles(simulation, polygon);
                    simulation.insert_no_flip(projected_point);
                    start_the_flips(simulation, polygon);
                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)){
                        /*Original insertion of Projection*/
                        custom_cdt.insert_no_flip(projected_point);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        cout<<"The projected_point steiner inserted 1st"<<endl;
                        break;
                    } 
                }
                /////////////////////////
            }
            
            if(count_obtuse_triangles(custom_cdt, polygon) == 0){
                CGAL::draw(custom_cdt);
                progress = false;
                break;
            }
        }
    }
    
    cout << "Mikelos has loops: " <<num<<" and we have "<<count_obtuse_triangles(custom_cdt, polygon)<<" obtuses"<<endl;
    cout<<"In is obtuse i have "<<obt<<" loops"<<endl;
    CGAL::draw(custom_cdt);
}

bool can_insert_on_projection(const Point_2& projected_point, const Polygon polygon) {
    //Check if you projected point is on side (boundary)
    if (is_point_inside_region(projected_point, polygon)) return true;
    else return false;
}

// Function to check if an edge is part of the boundary of the polygon
bool is_edge_in_boundary(const Point_2& p1, const Point_2& p2, const Polygon& polygon) {
    for (auto edge_it = polygon.edges_begin(); edge_it != polygon.edges_end(); ++edge_it) {
        if ((edge_it->source() == p1 && edge_it->target() == p2) ||
            (edge_it->source() == p2 && edge_it->target() == p1)) {
            return true;
        }
    }
    return false;
}

//Ιnsert Steiner points at circumcenters
void insert_circumcenter_centroid(Custom_CDT& custom_cdt, const Polygon& polygon) {
    
    cout<<"i am in  insert_steiner_points "<<endl;
    
    bool progress = true;
    while (progress ) {
        progress = false;
        
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            // Get the vertices of the current triangle
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();
            
            if(is_obtuse(p1, p2, p3)){ 
                //Compute the circumcenter of the triangle
                Face_handle triangleA = face;
                Point circumcenter = CGAL::circumcenter(p1, p2, p3);
                if (is_point_inside_region(circumcenter, polygon)){
                    
                    ////////////////////////////////
                    Custom_CDT simulation = custom_cdt;
                    if(is_convex(p1, p2, p3, circumcenter)){
                        cout<<"IS CONVEX"<<endl;

                        int initial_obtuse_count = count_obtuse_triangles(custom_cdt, polygon);
                        /*Simulate circumcenter insertion*/
                        simulation.insert_no_flip(circumcenter);
                        start_the_flips(simulation, polygon);
                        int final_obtuse_count = count_obtuse_triangles(simulation, polygon);

                        // Check if the flip resolved obtuse angles in the two faces
                        if (final_obtuse_count < initial_obtuse_count) {
                            /*Original circumcenter insertion*/
                            custom_cdt.insert_no_flip(circumcenter);
                            start_the_flips(custom_cdt, polygon);
                            progress = true;
                            final_obtuse_count = count_obtuse_triangles(custom_cdt, polygon);
                            //custom_cdt.flip(triangleA, i);
                            cout<<"Obtuse angles reduced: initial = "<<initial_obtuse_count<<", final = " <<final_obtuse_count<<endl;
                            break;
                        }
                    }
                    else{
                        Point_2 centroid = CGAL::centroid(p1, p2, p3);
                        //Custom_CDT simulation = custom_cdt;
                        bool insert_centroid = can_insert_centroid(custom_cdt, triangleA, centroid, polygon);
                        if(insert_centroid){
                            custom_cdt.insert(centroid);
                            start_the_flips(custom_cdt, polygon);
                            progress = true;
                            cout<<"The centroid steiner inserted"<<endl;
                            //CGAL::draw(custom_cdt);
                            break;
                        }
                    }
                } 
            }
        }
    }
}


//Just simulate if insert centroid (and auto flips) we reduce the obtuses
bool can_insert_centroid(Custom_CDT& custom_cdt, Face_handle& triangleA, const Point_2& centroid, const Polygon& polygon) {
    Custom_CDT simulation = custom_cdt;
    //Get the vertices of triangle A
    Point_2 p1 = triangleA->vertex(0)->point();
    Point_2 p2 = triangleA->vertex(1)->point();
    Point_2 p3 = triangleA->vertex(2)->point();

    int initial_obtuse_count = count_obtuse_triangles(custom_cdt, polygon);
    /*Simulate centroid insertion*/
    simulation.insert(centroid);
    start_the_flips(simulation, polygon);
    int final_obtuse_count = count_obtuse_triangles(simulation, polygon);

    //Check if the number of obtuse triangles decreased or stayed the same
    if (final_obtuse_count < initial_obtuse_count) {
        cout<<"Centroid Steiner insertion reduced obtuse angles: initial = "<<initial_obtuse_count<<", final = "<< final_obtuse_count<<endl;
        return true;
    } else {
        return false;
    }
}


bool is_point_inside_triangle(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& steiner_point) {
    //Determine the position of the steiner_point relative to the oriented circumcircle of the triangle
    CGAL::Oriented_side result = CGAL::side_of_oriented_circle(p1, p2, p3, steiner_point);

    //If the result is CGAL::ON_NEGATIVE_SIDE or CGAL::ON_ORIENTED_BOUNDARY, the point is inside the triangle or on the border
    return result == CGAL::ON_NEGATIVE_SIDE || result == CGAL::ON_ORIENTED_BOUNDARY;
    /*CGAL::ON_POSITIVE_SIDE means the point is outside the circumcircle (outside the triangle).
    CGAL::ON_NEGATIVE_SIDE means the point is inside the circumcircle.
    CGAL::ON_ORIENTED_BOUNDARY means the point lies on the circumcircle.*/
}

void start_the_flips(Custom_CDT& cdt, const Polygon& polygon){
    bool progress = true;
    while(progress){
        progress = false;
        for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
            CDT::Face_handle f1 = edge->first;
            int i = edge->second;
            CDT::Face_handle f2 = f1->neighbor(i);

            if (cdt.is_infinite(f1) || cdt.is_infinite(f2)) continue;
            
            Point_2 p1 = f1->vertex(cdt.ccw(i))->point(); // First vertex on the shared edge (Counter-Clock Wise)
            Point_2 p3 = f1->vertex(cdt.cw(i))->point();  // Second vertex on the shared edge (Clock Wise)
            Point_2 p2 = f1->vertex(i)->point();          // Opposite vertex in the first triangle
            //Mirror index gets the opposite vertex of the second triangle (f2)
            int mirror_index = cdt.mirror_index(f1, i);
            Point_2 p4 = f2->vertex(mirror_index)->point(); 
            //if the eddge is constraints or is on
            if (cdt.is_constrained(*edge) || is_edge_in_boundary(p1, p3, polygon)) continue;
            
            if(can_flip(p1,p2,p3,p4)){
                cdt.flip(f1, i);
                progress = true;
                cout<<"FLIPPED!"<<endl;
                break;
            }
            else continue;
        }
    }
}

//if 1 point is on the boundary
bool is_point_inside_region(const Point_2& point, const Polygon polygon) {

    //Check if the point is inside the polygon
    return (polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE) || (polygon.bounded_side(point) == CGAL::ON_BOUNDARY);
}

vector<Point_2> make_region_boundary(const vector<int>& region_bound_indx, const vector<Point_2>& points) {
    vector<Point_2> region_boundary;
    for (int index : region_bound_indx) {
        if (index >= 0 && index < points.size()) {
            region_boundary.push_back(points[index]);
        } else {
            //Handle error if the index is out of bounds
            std::cerr << "Error: Index " << index << " is out of bounds for the points vector." <<endl;
        }
    }
    return region_boundary;
}


Point_2 find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3) {
    // Calculate squared distances
    double ab2 = CGAL::to_double(squared_distance(v1, v2));
    double ac2 = CGAL::to_double(squared_distance(v1, v3));
    double bc2 = CGAL::to_double(squared_distance(v2, v3));

    // Check if the angle at v1 is obtuse
    if (ab2 + ac2 < bc2) {
        return v1; // obtuse angle at v1
    }
    // Check if the angle at v2 is obtuse
    if (ab2 + bc2 < ac2) {
        return v2; // obtuse angle at v2
    }
    // Check if the angle at v3 is obtuse
    if (ac2 + bc2 < ab2) {
        return v3; // obtuse angle at v3
    }

    throw std::logic_error("No obtuse angle found in the triangle.");
}

// Function to find the neighboring face that shares the edge (point1, point2) with the given face
Face_handle get_neighboring_face(Face_handle face, const Point_2& point1, const Point_2& point2) {
    // Loop through each edge of the given face to find the shared edge
    for (int i = 0; i < 3; ++i) {
        // Get the vertices of the current edge
        Point_2 v1 = face->vertex((i + 1) % 3)->point();
        Point_2 v2 = face->vertex((i + 2) % 3)->point();
        // Check if the edge matches the provided points
        if ((v1 == point1 && v2 == point2) || (v1 == point2 && v2 == point1)) {
            // Return the neighboring face that shares this edge
            return face->neighbor(i);
        }
    }

    // If no neighboring face is found, return a null handle
    return Face_handle();
}


Vertex_handle find_nearest_vertex(const Custom_CDT& cdt, const Point_2& query_point) {
    // Initialize with an invalid vertex handle
    Vertex_handle nearest_vertex = nullptr;
    double min_distance = std::numeric_limits<double>::max();

    for (auto v = cdt.finite_vertices_begin(); v != cdt.finite_vertices_end(); ++v) {
        double distance = CGAL::to_double(squared_distance(v->point(), query_point));
        if (distance < min_distance) {
            min_distance = distance;
            nearest_vertex = v;
        }
    }

    return nearest_vertex;
}

//if 2 points is on the boundary
bool are_points_on_polygon(const Polygon& polygon, const Point_2& p1, const Point_2& p2) {
    // Check if p1 and p2 are inside or on the boundary
    bool p1_on_polygon = (polygon.bounded_side(p1) == CGAL::ON_BOUNDED_SIDE ||
                          polygon.has_on_boundary(p1));
    bool p2_on_polygon = (polygon.bounded_side(p2) == CGAL::ON_BOUNDED_SIDE ||
                          polygon.has_on_boundary(p2));
    
    // Return true if both points are on the polygon (inside or on the boundary)
    return p1_on_polygon && p2_on_polygon;
}

CGAL::Segment_2<K> find_longest_edge(const Point_2& p1, const Point_2& p2, const Point_2& p3) {
    // Create segments for each edge of the triangle
    CGAL::Segment_2<K> edge1(p1, p2);
    CGAL::Segment_2<K> edge2(p2, p3);
    CGAL::Segment_2<K> edge3(p3, p1);

    // Compare the squared lengths of the segments to find the longest one
    auto length1 = CGAL::squared_distance(p1, p2);
    auto length2 = CGAL::squared_distance(p2, p3);
    auto length3 = CGAL::squared_distance(p3, p1);

    // Determine which edge is the longest
    if (length1 >= length2 && length1 >= length3) {
        return edge1;
    } else if (length2 >= length1 && length2 >= length3) {
        return edge2;
    } else {
        return edge3;
    }
}

    