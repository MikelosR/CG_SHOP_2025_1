#include "functions.h"
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"


using namespace boost::json;
using namespace std;
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

bool is_obtuse(const Point_2& a, const Point_2& b, const Point_2& c) {
    return (CGAL::angle(a, b, c) == CGAL::OBTUSE || 
            CGAL::angle(b, a, c) == CGAL::OBTUSE || 
            CGAL::angle(a, c, b) == CGAL::OBTUSE);
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

//DELETE
//Just check if an edge is part of the additional constrains
bool is_in_constraints(const std::pair<int, int>& edge, const vector<std::pair<int, int>>& constraints) {
    return std::find(constraints.begin(), constraints.end(), edge) != constraints.end() ||
           std::find(constraints.begin(), constraints.end(), std::make_pair(edge.second, edge.first)) != constraints.end();
           //checks for the edge in reverse order as well

}

//DELETE
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

                if (is_point_inside_region(midpoint, polygon)) {
                    int obtuses_before = count_obtuse_triangles(simulation, polygon);
                    simulation.insert_no_flip(midpoint);
                    start_the_flips(simulation, polygon);
                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)) {
                        custom_cdt.insert_no_flip(midpoint);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        break;
                    }
                }
            }
            if(count_obtuse_triangles(custom_cdt, polygon) == 0) progress = 0;
        }
    }
}


void insert_orthocenter(Custom_CDT& custom_cdt, const Polygon& polygon) {
    bool progress = true;
    while (progress) {
        progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();

            //Check if the triangle is obtuse
            if (is_obtuse(p1, p2, p3)) {
                Custom_CDT simulation = custom_cdt;
                Point_2 obtuse_vertex = find_obtuse_vertex(p1, p2, p3);

                //Calculate the orthocenter of the obtuse triangle
                Point_2 orthocenter = find_orthocenter(p1, p2, p3);

                //Verify if the orthocenter point can be inserted
                if (is_point_inside_region(orthocenter, polygon)) {
                    int obtuses_before = count_obtuse_triangles(simulation, polygon);
                    simulation.insert_no_flip(orthocenter);
                    start_the_flips(simulation, polygon);

                    //Check if insertion of orthocenter reduces obtuse triangles
                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)) {
                        custom_cdt.insert_no_flip(orthocenter);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        break;
                    }
                }
            }
            if(count_obtuse_triangles(custom_cdt, polygon) == 0) progress = 0;
        }
    }
}


void insert_projection(Custom_CDT& custom_cdt, const Polygon polygon){
    //Projection case
    //Progress flag to terminate the loop if after from 1 for loop without progress (no reduce the obtuses) to reminate the insert_projection
    bool progress = true;
    bool fill_boundary = true;
    while(progress){
        //If we end with the loops for the boundary faces, set the progress of main case false
        if(!fill_boundary) progress = false;
        for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face){
            Point_2 p1 = face->vertex(0)->point();
            Point_2 p2 = face->vertex(1)->point();
            Point_2 p3 = face->vertex(2)->point();
            Point_2 opposite1, opposite2;

            if(is_obtuse(p1, p2, p3)){
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
                //Copy the cdt for simulation
                Custom_CDT simulation = custom_cdt;
                bool insert_projection = is_point_inside_region(projected_point, polygon);
                
                int obtuses_before = count_obtuse_triangles(simulation, polygon);
                //Iterate all the faces of the boundary. If the point in side (not on) of the boundary is obtuse, insert projection
                if(fill_boundary && insert_projection){
                    if(is_face_on_boundary(custom_cdt, face) && !(polygon.bounded_side(obtuse_angle_vertex) == CGAL::ON_BOUNDARY)){
                        simulation.insert_no_flip(projected_point);
                        start_the_flips(simulation, polygon);
                        if (obtuses_before > count_obtuse_triangles(simulation, polygon)){
                            /*Original insertion of Projection*/
                            custom_cdt.insert_no_flip(projected_point);
                            start_the_flips(custom_cdt, polygon);
                        }
                    }
                }
                //main case, not only the boundary faces
                else if(insert_projection){
                    simulation.insert(projected_point);
                    start_the_flips(simulation, polygon);

                    if (obtuses_before > count_obtuse_triangles(simulation, polygon)){
                        /*Original insertion of Projection*/
                        custom_cdt.insert_no_flip(projected_point);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        break;
                    }
                    //Case the projection point reduce -1 the obtuse angle, but make new obtuse angle in the neighbor
                    //If this neighbor is boundary face (at least 2 points lies in the boundary polygon) and the projection point
                    //is part of this face but is in (not on) the boundary, add the projection
                    else if (obtuses_before == count_obtuse_triangles(simulation, polygon)) {
                    
                        for (auto sim_face = simulation.finite_faces_begin(); sim_face != simulation.finite_faces_end(); ++sim_face) {
                            if (simulation.is_infinite(sim_face)) continue;
                            
                            // Check if the projected point is part of a boundary face
                            if (sim_face->vertex(0)->point() == projected_point ||
                                sim_face->vertex(1)->point() == projected_point ||
                                sim_face->vertex(2)->point() == projected_point) {
                                if (is_face_on_boundary(simulation, sim_face)) {
                                    if(is_obtuse(sim_face->vertex(0)->point(), sim_face->vertex(1)->point(), sim_face->vertex(2)->point())){
                                        //Insert into the original triangulation
                                        custom_cdt.insert_no_flip(projected_point);
                                        start_the_flips(custom_cdt, polygon);
                                        progress = true;
                                        //"call" the fill boundary case
                                        fill_boundary = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            //when the case of iteration of boundary faces ends
            if((++face == custom_cdt.finite_faces_end()) && (fill_boundary)){
                fill_boundary = false;
            }
            --face;
            
            if(count_obtuse_triangles(custom_cdt, polygon) == 0){
                progress = false;
                break;
            }
            continue;
        }
    }
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
                    Custom_CDT simulation = custom_cdt;
                    if(is_convex(p1, p2, p3, circumcenter)){
                        int initial_obtuse_count = count_obtuse_triangles(custom_cdt, polygon);
                        /*Simulate circumcenter insertion*/
                        simulation.insert_no_flip(circumcenter);
                        start_the_flips(simulation, polygon);
                        int final_obtuse_count = count_obtuse_triangles(simulation, polygon);

                        //Check if the flip resolved obtuse angles in the two faces
                        if (final_obtuse_count < initial_obtuse_count) {
                            /*Original circumcenter insertion*/
                            custom_cdt.insert_no_flip(circumcenter);
                            start_the_flips(custom_cdt, polygon);
                            progress = true;
                            final_obtuse_count = count_obtuse_triangles(custom_cdt, polygon);
                            break;
                        }
                    }
                }
                else{
                    //Cetroid case (if circumcenter went out of bounds)
                    Point_2 centroid = CGAL::centroid(p1, p2, p3);
                    bool insert_centroid = can_insert_centroid(custom_cdt, triangleA, centroid, polygon);
                    //Just check if the insertion of centroid has a benefit
                    if(insert_centroid){
                        custom_cdt.insert(centroid);
                        start_the_flips(custom_cdt, polygon);
                        progress = true;
                        break;
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
    if (final_obtuse_count < initial_obtuse_count) return true;
    else return false;
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
            //if the edge is constraints or is on boundary
            if (cdt.is_constrained(*edge) || is_edge_in_boundary(p1, p3, polygon)) continue;
            
            if(can_flip(p1, p2, p3, p4)){
                cdt.flip(f1, i);
                progress = true;
                break;
            }
            else continue;
        }
    }
}

//If 1 point is on the boundary
bool is_point_inside_region(const Point_2& point, const Polygon polygon) {

    //Check if the point is inside the polygon
    return (polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE) || (polygon.bounded_side(point) == CGAL::ON_BOUNDARY);
}


Point_2 find_obtuse_vertex(const Point_2& v1, const Point_2& v2, const Point_2& v3) {
    //Calculate squared distances
    double ab2 = CGAL::to_double(squared_distance(v1, v2));
    double ac2 = CGAL::to_double(squared_distance(v1, v3));
    double bc2 = CGAL::to_double(squared_distance(v2, v3));

    if (ab2 + ac2 < bc2) return v1; // obtuse angle at v1
    if (ab2 + bc2 < ac2) return v2; // obtuse angle at v2
    if (ac2 + bc2 < ab2) return v3; // obtuse angle at v3

    throw std::logic_error("No obtuse angle found in the triangle.");
    //return Point_2(0, 0);
}

CGAL::Segment_2<K> find_longest_edge(const Point_2& p1, const Point_2& p2, const Point_2& p3) {
    //Create segments for each edge of the triangle
    CGAL::Segment_2<K> edge1(p1, p2);
    CGAL::Segment_2<K> edge2(p2, p3);
    CGAL::Segment_2<K> edge3(p3, p1);

    //Compare the squared lengths of the segments to find the longest one
    auto length1 = CGAL::squared_distance(p1, p2);
    auto length2 = CGAL::squared_distance(p2, p3);
    auto length3 = CGAL::squared_distance(p3, p1);

    //Determine which edge is the longest
    if (length1 >= length2 && length1 >= length3) {
        return edge1;
    } else if (length2 >= length1 && length2 >= length3) {
        return edge2;
    } else {
        return edge3;
    }
}

Point_2 find_medial_of_longest_side(const Point_2& p1, const Point_2& p2, const Point_2& p3) {
    // Calculate the squared lengths of each side
    double length1_sq = CGAL::to_double(CGAL::squared_distance(p1, p2)); // Side p1 - p2
    double length2_sq = CGAL::to_double(CGAL::squared_distance(p2, p3)); // Side p2 - p3
    double length3_sq = CGAL::to_double(CGAL::squared_distance(p1, p3)); // Side p1 - p3

    // Find the longest side
    if (length1_sq >= length2_sq && length1_sq >= length3_sq) {
        // Longest side is p1 - p2
        return CGAL::midpoint(p1, p2);
    } else if (length2_sq >= length1_sq && length2_sq >= length3_sq) {
        // Longest side is p2 - p3
        return CGAL::midpoint(p2, p3);
    } else {
        // Longest side is p1 - p3
        return CGAL::midpoint(p1, p3);
    }
}

Point_2 find_orthocenter(const Point_2& p1, const Point_2& p2, const Point_2& p3) {
    //Calculate slopes and midpoints of sides
    Line_2 line1(p1, p2);
    Line_2 line2(p2, p3);
    Line_2 line3(p1, p3);

    //Find the orthogonal bisectors
    Line_2 perp1 = line1.perpendicular(p3);
    Line_2 perp2 = line2.perpendicular(p1);

    //Intersection of two perpendiculars gives the orthocenter
    CGAL::Object result = CGAL::intersection(perp1, perp2);
    const Point_2* orthocenter = CGAL::object_cast<Point_2>(&result);
    
    //Ensure result validity
    if (orthocenter) return *orthocenter;
    else throw std::runtime_error("Orthocenter calculation failed.");
}

//Found if this face is on the boundary!!
bool is_face_on_boundary(const Custom_CDT& cdt, Face_handle face) {
    //Loop through each edge of the face
    for (int i = 0; i < 3; ++i) {
        Custom_CDT::Edge edge(face, i);
        if (cdt.is_infinite(face->neighbor(i))) {
            return true; //Edge is on the boundary
        }
    }
    return false;
}

//Function to check if a vertex in custom_cdt is a Steiner point
bool is_steiner_point(Vertex_handle vertex, const vector<Point_2>& original_points) {
    const Point_2& v_point = vertex->point();
    //Check if v_point is in the original_points vector
    return std::none_of(original_points.begin(), original_points.end(), [&v_point](const Point_2& p) { return p == v_point; });
}


void output(value jv, Custom_CDT custom_cdt, vector<CGAL::Point_2<CGAL::Epeck>> points) {
    // Get instance_uid from input JSON
    std::string instance_uid = jv.as_object().at("instance_uid").as_string().c_str();

    // Prepare steiner points lists
    vector<std::string> steiner_points_x;
    vector<std::string> steiner_points_y;

    for (auto vertex = custom_cdt.finite_vertices_begin(); vertex != custom_cdt.finite_vertices_end(); ++vertex) {
        // if a vertex is not found in the initial points then it must be a steiner point
        if (is_steiner_point(vertex, points)) {  
            const Point_2& p = vertex->point();
            steiner_points_x.push_back(convert_to_string(p.x()));  // Store x-coordinate as string
            steiner_points_y.push_back(convert_to_string(p.y()));  // Store y-coordinate as string
        }
    }

    // Convert steiner_points_x and steiner_points_y to JSON-compatible format
    boost::json::array steiner_points_x_json;
    boost::json::array steiner_points_y_json;

    for (const auto& x : steiner_points_x) {
        steiner_points_x_json.push_back(boost::json::value(x));
    }

    for (const auto& y : steiner_points_y) {
        steiner_points_y_json.push_back(boost::json::value(y));
    }

    // Prepare edges list
    vector<std::pair<int, int>> edges;
    std::map<Vertex_handle, int> vertex_index_map;
    int index = 0;

    // Map vertices to unique indices for edge representation
    for (auto vertex = custom_cdt.finite_vertices_begin(); vertex != custom_cdt.finite_vertices_end(); ++vertex) {
        vertex_index_map[vertex] = index++;
    }

    for (auto edge = custom_cdt.finite_edges_begin(); edge != custom_cdt.finite_edges_end(); ++edge) {
        Vertex_handle v1 = edge->first->vertex((edge->second + 1) % 3);
        Vertex_handle v2 = edge->first->vertex((edge->second + 2) % 3);
        edges.emplace_back(vertex_index_map[v1], vertex_index_map[v2]);
    }

    // Convert edges to json-compatible format
    boost::json::array edges_json;
    for (const auto& edge : edges) {
        boost::json::array edge_pair = {edge.first, edge.second};
        edges_json.push_back(edge_pair);
    }

    // make the json file look pretty (not in a single line)
    std::ostringstream json_output;
    json_output << "{\n";
    json_output << "  \"content_type\": \"CG_SHOP_2025_Solution\",\n";
    json_output << "  \"instance_uid\": \"" << instance_uid << "\",\n";
    json_output << "  \"steiner_points_x\": " << boost::json::serialize(steiner_points_x_json) << ",\n";
    json_output << "  \"steiner_points_y\": " << boost::json::serialize(steiner_points_y_json) << ",\n";
    json_output << "  \"edges\": " << boost::json::serialize(edges_json) << "\n";
    json_output << "}\n";

    // write json to file
    std::ofstream output_file("solution_output.json");
    output_file << json_output.str();
    output_file.close();

    cout << "Solution JSON file written as solution_output.json" << endl;
    return;
}

//Convert coord into string like franction. If it is int, return string like int
std::string convert_to_string(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    std::ostringstream oss;
     // Check if the coordinate is an integer
    if (exact_coord.get_den() == 1) {
        oss << coord;
        return oss.str();   
    }
    oss << exact_coord.get_num() << "/" << exact_coord.get_den();
    
    // Return the formatted string
    return oss.str();
}


//Alternative is_steiner_point
/*bool is_steiner_point(const Custom_CDT::Vertex_handle& vertex, const std::vector<Point_2>& initial_points) {
    return std::find(initial_points.begin(), initial_points.end(), vertex->point()) == initial_points.end();
}*/

