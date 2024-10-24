#include "functions.h"
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"
//#include <CGAL/Polygon_2_algorithms.h>

using namespace boost::json;
//using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
typedef CGAL::Polygon_2<K> Polygon;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Point_2 = K::Point_2;
using Line_2 = K::Line_2;

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
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    std::string json_str((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    jv = parse(json_str);
}

//Just check if an edge is part of the additional constrains
bool is_in_constraints(const std::pair<int, int>& edge, const std::vector<std::pair<int, int>>& constraints) {
    return std::find(constraints.begin(), constraints.end(), edge) != constraints.end() ||
           std::find(constraints.begin(), constraints.end(), std::make_pair(edge.second, edge.first)) != constraints.end();
           //checks for the edge in reverse order as well

}

//Just check if an edge is part of the region boundary
bool is_in_region_boundary(const std::pair<int, int>& edge, const std::vector<int>& boundary) {
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
int count_obtuse_triangles(CDT& cdt) {
    int obtuse_count = 0;
    for (auto fit = cdt.finite_faces_begin(); fit != cdt.finite_faces_end(); ++fit) {
        Point_2 p1 = fit->vertex(0)->point();
        Point_2 p2 = fit->vertex(1)->point();
        Point_2 p3 = fit->vertex(2)->point();
        if (is_obtuse(p1, p2, p3)) {
            obtuse_count++;
        }
    }
    return obtuse_count;
}

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point_2& p1, const Point_2& p2, const Point_2& p3, const Point_2& p4) {
    return CGAL::orientation(p2, p1, p4) == CGAL::LEFT_TURN &&
           CGAL::orientation(p4, p3, p2) == CGAL::LEFT_TURN;
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

/*void insert_points_within_boundary(CDT& cdt, const std::vector<Point_2>& points, const std::vector<Point_2>& boundary_points) {
    // Construct the polygon from the boundary points
    Polygon boundary_polygon(boundary_points.begin(), boundary_points.end());

    // Iterate through all the points and insert only those that lie inside the boundary
    for (const auto& point : points) {
        if (boundary_polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE) {
            cdt.insert(point); // Insert point if it is inside the polygon
        }
    }

    // Insert the boundary edges as constraints
    for (size_t i = 0; i < boundary_points.size(); ++i) {
        const Point_2& p1 = boundary_points[i];
        const Point_2& p2 = boundary_points[(i + 1) % boundary_points.size()];
        cdt.insert_constraint(p1, p2); // Insert each edge as a constraint
    }
}*/

//Ιnsert Steiner points at circumcenters
void insert_steiner_points(Custom_CDT& custom_cdt, const std::vector<Point_2>& points, const std::vector<int>& region_boundary, const Polygon polygon) {
    //Iterate over all faces in the triangulation
    std::cout<<"i am in  insert_steiner_points "<<std::endl;
    bool progress = true;
    int iteration_count = 0;            // Counter to avoid infinite loop
    const int max_iterations = 1000;    // Set a reasonable maximum iteration limit

    while (progress && iteration_count < max_iterations) {
        progress = false;
        iteration_count++;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            // Get the vertices of the current triangle
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            if(is_obtuse(p1, p2, p3)){ 
                //Compute the circumcenter of the triangle
                Custom_CDT simulation = custom_cdt;
                Custom_CDT::Face_handle triangleA = face;
                Point circumcenter = CGAL::circumcenter(p1, p2, p3);
                if (is_point_inside_region(circumcenter, polygon) && 
                    !is_obtuse(p1, p2, circumcenter) && 
                    !is_obtuse(p2, p3, circumcenter) && 
                    !is_obtuse(p3, p1, circumcenter)){

                    bool insert_circumcenter = can_insert_circumcenter(simulation, triangleA, circumcenter);
                    if (insert_circumcenter) {
                        custom_cdt.insert_no_flip(circumcenter);
                        progress = true;
                        std::cout << "The circumcenter Steiner inserted" << std::endl;
                        //CGAL::draw(custom_cdt);
                        break;
                    }
                    else{
                        Point_2 centroid = CGAL::centroid(p1, p2, p3);
                        bool insert_centroid = can_insert_centroid(simulation, triangleA, centroid);
                        if(insert_centroid){
                            custom_cdt.insert_no_flip(centroid);
                            progress = true;
                            std::cout<<"The centroid steiner inserted"<<std::endl;
                            //CGAL::draw(custom_cdt);
                            break;
                        }
                    }
                }
                //Projection case
                Point_2 obtuse_angle_vertex = find_obtuse_vertex(p1, p2, p3);
                Point_2 opposite1, opposite2;

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
                // Create a line from the edge opposite the obtuse angle
                Line_2 line(opposite1, opposite2);
                // Project the obtuse angle vertex onto the line
                Point_2 projected_point = line.projection(obtuse_angle_vertex);
                bool insert_projection = can_insert_on_projection(projected_point, polygon);

                // Find the neighboring triangle (Triangle B) that shares the edge with Triangle A
                Custom_CDT::Face_handle neighbor = get_neighboring_face(face, opposite1, opposite2);
                
                if (neighbor != nullptr) {
                    //Get the fourth point from the neighboring triangle (Triangle B)
                    Point_2 p4;
                    for (int i = 0; i < 3; ++i) {
                        //maybe it wrong this way....
                        if (neighbor->vertex(i)->point() != opposite1 && neighbor->vertex(i)->point() != opposite2) {
                            p4 = neighbor->vertex(i)->point();
                            break;
                        }
                    }
                    
                    if(insert_projection){
                        custom_cdt.insert_no_flip(projected_point);
                        progress = true;
                        std::cout<<"The projected_point steiner inserted"<<std::endl;
                        //CGAL::draw(custom_cdt);
                        // Get the vertex handles for the projected point and p4
                        Custom_CDT::Vertex_handle vh_projected = custom_cdt.nearest_vertex(projected_point);
                        Custom_CDT::Vertex_handle vh_p4 = custom_cdt.nearest_vertex(p4);
                        //If the edge already exists, no need to add
                        if (custom_cdt.is_edge(vh_projected, vh_p4)) {
                            std::cout << "Edge already exists between projected point and p4." << std::endl;
                        } else {
                            //Add edge [projected_point, p4] to maintain triangulation
                            custom_cdt.insert_constraint(projected_point, p4);
                            std::cout << "Inserted edge between projected point and p4." << std::endl;
                        }
                        break;
                    }
                }
            }
        }
    }
}

bool can_insert_on_projection(const Point_2& projected_point, const Polygon polygon) {
    //Check if you projected point is on side (boundary)
    if (is_point_inside_region(projected_point, polygon)) return true;
    else return false;
}
                 

bool can_insert_circumcenter(Custom_CDT& custom_cdt, Custom_CDT::Face_handle& triangleA, const Point_2& circumcenter) {
    // Count the number of obtuse triangles before inserting the Steiner point
    int initial_obtuse_count = count_obtuse_triangles(custom_cdt);
    
    // Insert the circumcenter without flipping any edges
    custom_cdt.insert_no_flip(circumcenter);

    //For every edge of the triangleA, find this edge with the neighbor has the circumcenter as point
    for (int i = 0; i < 3; ++i) {
        // Find the neighboring triangle B
        Custom_CDT::Face_handle triangleB = triangleA->neighbor(i);

        // If the neighboring triangle is not infinite and contains the Steiner point, we have found triangle B
        if (!custom_cdt.is_infinite(triangleB)) {
            // Check if the neighboring triangle B contains the circumcenter as a vertex
            if (triangleB->vertex(0)->point() == circumcenter ||
                triangleB->vertex(1)->point() == circumcenter ||
                triangleB->vertex(2)->point() == circumcenter) {
                
                //you are ensuring that p1 and p2 correspond to the vertices that are shared between triangleA and triangleB,
                //and p3 is the third vertex of triangleA.
                //triangleA->index(triangleB); returns the local index (0, 1, or 2) of the vertex in triangleA opposite to the 
                //shared edge between triangleA and triangleB
                //triangleA->index(triangleB) gives you the index of the vertex in triangleA that is not part of the shared edge.
                int indexA1 = triangleA->index(triangleB);
                //the vertices that form the shared edge (indexA2, indexA3)
                int indexA2 = (indexA1 + 1) % 3;
                int indexA3 = (indexA1 + 2) % 3;

                // Vertices of Triangle A
                Point_2 p1 = triangleA->vertex(indexA1)->point();
                Point_2 p2 = triangleA->vertex(indexA2)->point();
                Point_2 p3 = triangleA->vertex(indexA3)->point();
                // Vertex of Triangle B not in Triangle A (which is the circumcenter)
                Point_2 p4 = circumcenter;

                // Check if flipping the edge can reduce the number of obtuse angles
                if (can_flip(p1, p2, p3, p4)) {
                    custom_cdt.flip(triangleA, i); // Perform the flip on the shared edge (i)
                    
                    // Count the number of obtuse triangles after the flip
                    int final_obtuse_count = count_obtuse_triangles(custom_cdt);
                    
                    // Check if the flip resolved obtuse angles in the two faces
                    if (final_obtuse_count < initial_obtuse_count) {
                        std::cout << "Obtuse angles reduced: initial = " << initial_obtuse_count
                                  << ", final = " << final_obtuse_count << std::endl;
                        return true; // Successful insertion with obtuse angles resolved
                    }
                }
                break;
            }
        }
    }
    return false; // Flip wasn't possible or didn't resolve obtuse angles
}


bool can_insert_centroid(Custom_CDT& custom_cdt, Custom_CDT::Face_handle& triangleA, const Point_2& centroid) {
    // Get the vertices of triangle A
    Point_2 p1 = triangleA->vertex(0)->point();
    Point_2 p2 = triangleA->vertex(1)->point();
    Point_2 p3 = triangleA->vertex(2)->point();

    // Check the number of obtuse triangles before inserting the centroid
    int initial_obtuse_count = 0;
    if (is_obtuse(p1, p2, p3)) {
        initial_obtuse_count = 1;
    }

    // Insert the centroid without flipping
    custom_cdt.insert_no_flip(centroid);

    // Now we have three new triangles: (p1, p2, centroid), (p2, p3, centroid), and (p3, p1, centroid)
    int final_obtuse_count = 0;
    if (is_obtuse(p1, p2, centroid)) final_obtuse_count++;
    if (is_obtuse(p2, p3, centroid)) final_obtuse_count++;
    if (is_obtuse(p3, p1, centroid)) final_obtuse_count++;

    // Check if the number of obtuse triangles decreased or stayed the same
    if (final_obtuse_count < initial_obtuse_count) {
        std::cout << "Centroid Steiner insertion reduced obtuse angles: initial = "
                  << initial_obtuse_count << ", final = " << final_obtuse_count << std::endl;
        return true; // Successful insertion with fewer obtuse angles
    } else {
        std::cout << "Centroid Steiner insertion did not reduce obtuse angles: initial = "
                  << initial_obtuse_count << ", final = " << final_obtuse_count << std::endl;
        return false; // Insertion did not reduce obtuse angles
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

void start_the_flips(Custom_CDT& cdt, const std::vector<Point_2>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary){
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        CDT::Face_handle f1 = edge->first;
        int i = edge->second;
        CDT::Face_handle f2 = f1->neighbor(i);

        if (cdt.is_infinite(f1) || cdt.is_infinite(f2)) {
            continue;
        }

        Point_2 p1 = f1->vertex(cdt.ccw(i))->point(); // First vertex on the shared edge (Counter-Clock Wise)
        Point_2 p3 = f1->vertex(cdt.cw(i))->point();  // Second vertex on the shared edge (Clock Wise)
        Point_2 p2 = f1->vertex(i)->point();          // Opposite vertex in the first triangle
        //Mirror index gets the opposite vertex of the second triangle (f2)
        int mirror_index = cdt.mirror_index(f1, i);
        Point_2 p4 = f2->vertex(mirror_index)->point(); 

        // Check if edge belongs to additional constraints or region boundary
        /*std::find(points.begin(), points.end(), p1): This function searches through the points vector to find the position of p1
        If the value is found in the range, returns an iterator to its position.
        if the value is not found, returns an iterator to the last position.*/
        //std::distance(InputIterator first, InputIterator last): This calculates the distance (number of steps) from the beginning of the points vector 
        //to the position where p1 is found. This gives us the index of p1 in the points vector
        
        int idx1 = std::distance(points.begin(), std::find(points.begin(), points.end(), p1));
        int idx2 = std::distance(points.begin(), std::find(points.begin(), points.end(), p3));

        if (is_in_constraints({idx1, idx2}, additional_constraints) || 
            is_in_region_boundary({idx1, idx2}, region_boundary)) {
            continue;
        }

        if(can_flip(p1,p2,p3,p4)){
            cdt.flip(f1, i);
            std::cout << "FLIPPED!" << std::endl;
            //CGAL::draw(cdt);
        }
        else{ 
            //CGAL::draw(cdt);
            continue; 
        }
    }
}

bool is_point_inside_region(const Point_2& point, const Polygon polygon) {

    //Check if the point is inside the polygon
    return (polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE) || (polygon.bounded_side(point) == CGAL::ON_BOUNDARY);
}

std::vector<Point_2> make_region_boundary(const std::vector<int>& region_bound_indx, const std::vector<Point_2>& points) {
    std::vector<Point_2> region_boundary;
    for (int index : region_bound_indx) {
        if (index >= 0 && index < points.size()) {
            region_boundary.push_back(points[index]);
        } else {
            //Handle error if the index is out of bounds
            std::cerr << "Error: Index " << index << " is out of bounds for the points vector." << std::endl;
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
Custom_CDT::Face_handle get_neighboring_face(Custom_CDT::Face_handle face, const Point_2& point1, const Point_2& point2) {
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
    return Custom_CDT::Face_handle();
}


