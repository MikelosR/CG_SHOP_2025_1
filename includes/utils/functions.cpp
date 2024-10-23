#include "functions.h"
#include "includes/utils/Custom_Constrained_Delaunay_triangulation_2.h"
//#include <CGAL/Polygon_2_algorithms.h>

using namespace boost::json;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
typedef CGAL::Polygon_2<K> Polygon;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;


//Calculate squared distance between two points
double squared_distance(const Point& a, const Point& b) {
    return CGAL::squared_distance(a, b);
}

//Check if a triangle is obtuse
bool is_obtuse(const Point& a, const Point& b, const Point& c) {
    double ab2 = squared_distance(a, b);
    double ac2 = squared_distance(a, c);
    double bc2 = squared_distance(b, c);
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
    for (size_t i = 0; i < boundary.size(); ++i) {
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
        Point p1 = fit->vertex(0)->point();
        Point p2 = fit->vertex(1)->point();
        Point p3 = fit->vertex(2)->point();
        if (is_obtuse(p1, p2, p3)) {
            obtuse_count++;
        }
    }
    return obtuse_count;
}

//Return true if 2 faces (two triangles) form a convex polygon
bool is_convex(const Point& p1, const Point& p2, const Point& p3, const Point& p4) {
    return CGAL::orientation(p2, p1, p4) == CGAL::LEFT_TURN &&
           CGAL::orientation(p4, p3, p2) == CGAL::LEFT_TURN;
}

//Return true if approves the flip
bool can_flip(const Point& p1, const Point& p2, const Point& p3, const Point& p4){
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

void insert_points_within_boundary(CDT& cdt, const std::vector<Point>& points, const std::vector<Point>& boundary_points) {
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
        const Point& p1 = boundary_points[i];
        const Point& p2 = boundary_points[(i + 1) % boundary_points.size()];
        cdt.insert_constraint(p1, p2); // Insert each edge as a constraint
    }
}

//Ιnsert Steiner points at circumcenters
void insert_steiner_points(Custom_CDT& custom_cdt, const std::vector<Point>& points, const std::vector<int>& region_boundary) {
    //Iterate over all faces in the triangulation
    std::cout<<"i am in  insert_steiner_points "<<std::endl;
    for (auto face = custom_cdt.finite_faces_begin(); face != custom_cdt.finite_faces_end(); ++face) {
        // Get the vertices of the current triangle
        Point p1 = face->vertex(0)->point();
        Point p2 = face->vertex(1)->point();
        Point p3 = face->vertex(2)->point();

        if(is_obtuse(p1,p2,p3)){ 
            //Compute the circumcenter of the triangle
            Point circumcenter = CGAL::circumcenter(p1, p2, p3);
            bool circmcentr = circumcenter_steiner(p1,p2,p3,circumcenter, points, region_boundary);
            if(circmcentr){
                custom_cdt.insert_no_flip(circumcenter);
                std::cout<<"The circumcenter steiner inserted"<<std::endl;
                CGAL::draw(custom_cdt);
            }
            else{
                //std::cout<<"The circumcenter steiner is outside of the triangle"<<std::endl;
                //Calculate the centroid
                Point centroid = CGAL::centroid(p1, p2, p3);
                bool cntroid = centroid_steiner(p1,p2,p3,centroid);
                if(cntroid){
                    custom_cdt.insert_no_flip(centroid);
                    std::cout<<"The centroid steiner inserted"<<std::endl;
                    CGAL::draw(custom_cdt);
                }
            }
        }
    }
}

bool centroid_steiner(const Point& p1, const Point& p2, const Point& p3, const Point& centroid){
    std::cout<<"The centroid CHECKED"<<std::endl;
    if (!is_obtuse(p1, p2, centroid) && !is_obtuse(p2, p3, centroid) && !is_obtuse(p3, p1, centroid)) {
        return true;
    }
    return false;
    
}   

bool circumcenter_steiner(const Point& p1, const Point& p2, const Point& p3, const Point& circumcenter, const std::vector<Point>& points,  const std::vector<int>& region_boundary ){
    
    //Check if the circumcenter is within the region boundary before inserting
    if (is_point_inside_region(circumcenter, region_boundary, points) && !is_obtuse(p1, p2, circumcenter) && 
        !is_obtuse(p2, p3, circumcenter) && !is_obtuse(p3, p1, circumcenter)) {
        
        //Check if the circumcenter (Steiner point) lies inside or on the boundary of the triangle
        if (is_point_inside_triangle(p1, p2, p3, circumcenter)) {
            std::cout<<"i am trying to insert the steiner (circumcenter): "<<circumcenter<<std::endl;           
            return true;
        }
    }
    return false;  
}

bool is_point_inside_triangle(const Point& p1, const Point& p2, const Point& p3, const Point& steiner_point) {
    //Determine the position of the steiner_point relative to the oriented circumcircle of the triangle
    CGAL::Oriented_side result = CGAL::side_of_oriented_circle(p1, p2, p3, steiner_point);

    //If the result is CGAL::ON_NEGATIVE_SIDE or CGAL::ON_ORIENTED_BOUNDARY, the point is inside the triangle or on the border
    return result == CGAL::ON_NEGATIVE_SIDE || result == CGAL::ON_ORIENTED_BOUNDARY;
    /*CGAL::ON_POSITIVE_SIDE means the point is outside the circumcircle (outside the triangle).
    CGAL::ON_NEGATIVE_SIDE means the point is inside the circumcircle.
    CGAL::ON_ORIENTED_BOUNDARY means the point lies on the circumcircle.*/
}

void start_the_flips(Custom_CDT& cdt, const std::vector<Point>& points, const std::vector<std::pair<int, int>>& additional_constraints, const std::vector<int>& region_boundary){
    for (auto edge = cdt.finite_edges_begin(); edge != cdt.finite_edges_end(); ++edge) {
        CDT::Face_handle f1 = edge->first;
        int i = edge->second;
        CDT::Face_handle f2 = f1->neighbor(i);

        if (cdt.is_infinite(f1) || cdt.is_infinite(f2)) {
            continue;
        }

        Point p1 = f1->vertex(cdt.ccw(i))->point(); // First vertex on the shared edge (Counter-Clock Wise)
        Point p3 = f1->vertex(cdt.cw(i))->point();  // Second vertex on the shared edge (Clock Wise)
        Point p2 = f1->vertex(i)->point();          // Opposite vertex in the first triangle
        //Mirror index gets the opposite vertex of the second triangle (f2)
        int mirror_index = cdt.mirror_index(f1, i);
        Point p4 = f2->vertex(mirror_index)->point(); 

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

bool is_point_inside_region(const Point& point, const std::vector<int>& region_boundary, const std::vector<Point>& points) {
    // Create a Polygon from the points in the region boundary
    Polygon polygon;
    for (int index : region_boundary) {
        polygon.push_back(points[index]);
    }

    //Check if the point is inside the polygon
    return polygon.bounded_side(point) == CGAL::ON_BOUNDED_SIDE ||
           polygon.bounded_side(point) == CGAL::ON_BOUNDARY;
}

std::vector<Point> make_region_boundary(const std::vector<int>& region_bound_indx, const std::vector<Point>& points) {
    std::vector<Point> region_boundary;
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