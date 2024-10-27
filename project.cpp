#include "includes/utils/functions.h"

using namespace boost::json;
using namespace std;
using K = CGAL::Exact_predicates_exact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;
using Point_2 = K::Point_2;
using Vertex_handle = CDT::Vertex_handle;


int main() {
    value jv;
    read_json("tests/pdf.json", jv);                  //5 obtuses 76%  success
    //read_json("tests/instance_test_1.json", jv);      //0 obtuses 100% success
    //read_json("tests/instance_test_2.json", jv);      //0 obtuses 100% success
    //read_json("tests/instance_test_3.json", jv);      //0 obtuse  100% success
    //read_json("tests/instance_test_4.json", jv);      //0 obtuses 100% success
    //read_json("tests/instance_test_5.json", jv);      //0 obtuses 100% success
    //read_json("tests/instance_test_6.json", jv);      //1 obtuse   75% success
    //read_json("tests/instance_test_7.json", jv);      //1 obtuse   75% success
    //read_json("tests/instance_test_8.json", jv);      //1 obtuse   80% success
    //read_json("tests/instance_test_9.json", jv);      //1 obtuse   86% success
    //read_json("tests/instance_test_10.json", jv);     //1 obtuse   66% success
    //read_json("tests/instance_test_11.json", jv);     //6 obtuses  40% success
    //read_json("tests/instance_test_12.json", jv);     //0 obtuses 100% success
    //read_json("tests/instance_test_13.json", jv);     //4 obtuses  33% success
    //read_json("tests/instance_test_14.json", jv);     //2 obtuses  60% success
    //read_json("tests/instance_test_15.json", jv);     //2 obtuses  60% success
    //read_json("tests/instance_test_16.json", jv);     //1 obtuse   80% success

    //read_json("tests/instance_1_1.json", jv);         //0 obtuses 100% success
    //read_json("tests/instance_2_1.json", jv);         //0 obtuses 100% success
    //read_json("tests/instance_3_1.json", jv);         //0 obtuses 100% success
    //read_json("tests/instance_4_1.json", jv);         //0 obtuses 60% success
    //read_json("tests/instance_5_1.json", jv);         //0 obtuses 10% success    
    //read_json("tests/instance_6_1.json", jv);         //0 obtuses 66% success
    //read_json("tests/instance_7_1.json", jv);         //0 obtuses 83% success
     
    vector<Point_2> points;
    vector<std::pair<int, int>> additional_constraints;
    vector<int> region_boundary;
    Polygon polygon;

    if (jv.is_object()) {
        const auto& obj = jv.as_object();
        const auto& x_array = obj.at("points_x").as_array();
        const auto& y_array = obj.at("points_y").as_array();
        const auto& boundary_array = obj.at("region_boundary").as_array();
        const auto& constraints_array = obj.at("additional_constraints").as_array();

        for (int i = 0; i < x_array.size(); ++i) {
            double x = x_array[i].is_double() ? x_array[i].as_double() : static_cast<double>(x_array[i].as_int64());
            double y = y_array[i].is_double() ? y_array[i].as_double() : static_cast<double>(y_array[i].as_int64());
            points.emplace_back(x, y);
        }

        for (const auto& idx : boundary_array) {
            region_boundary.push_back(idx.as_int64());
        }

        //Add the additional constraints in vector
        for (const auto& constraint : constraints_array) {
            int idx1 = constraint.as_array()[0].as_int64();
            int idx2 = constraint.as_array()[1].as_int64();
            if (idx1 < points.size() && idx2 < points.size()) {
                additional_constraints.emplace_back(idx1, idx2);
            }
        }
    }

    //Create a polygon from region boundary
    for (int index : region_boundary) {
        polygon.push_back(points[index]);
    }

    //Make the cdt
    Custom_CDT custom_cdt;
    vector<Vertex_handle> vertex_handles;
    for (const auto& point : points) {
        vertex_handles.push_back(custom_cdt.insert(point));
    }

    //Insert additional constraints
    for (const auto& constraint : additional_constraints) {
        custom_cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }
    
    //If we have progress (reduce obtuses) run again
    bool progress = true;
    while(progress){
        progress = false;
        int start = count_obtuse_triangles(custom_cdt, polygon);
        int num_obtuses_before = start;
        cout<<"Number of obtuse triangles before FLIPS: " <<num_obtuses_before<<endl;
        /*Flips*/
        start_the_flips(custom_cdt, polygon);
        int num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout <<"Number of obtuse triangles after FLIPS : "<<num_obtuses_after<<endl;
        CGAL::draw(custom_cdt);
        if(num_obtuses_after == 0) break;

        //Circumcenter - Centroid
        num_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles before inserting CIRCUMCENTER and CENTROIDS: "<<num_obtuses_before<<endl;
        /*Steiner Points circumcenter and centroids*/
        insert_circumcenter_centroid(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout << "Number of obtuse triangles after inserting CIRCUMCENTER and CENTROIDS: "<<num_obtuses_after<<endl;
        //CGAL::draw(custom_cdt);
        if(num_obtuses_after == 0) break;

        //Midpoint
        num_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles inserting MIDPOINT: " <<num_obtuses_before<<endl;
        /*Steiner Points Midpoint*/
        insert_midpoint(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles after inserting MIDPOINT: "<<num_obtuses_after<<endl;
        //CGAL::draw(custom_cdt);
        if(num_obtuses_after == 0) break;

        //Projection
        num_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles inserting PROJECTION: " <<num_obtuses_before<<endl;
        /*Steiner Points Projections*/
        insert_projection(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles after inserting PROJECTION: "<<num_obtuses_after<<endl;
        //CGAL::draw(custom_cdt);
        if(num_obtuses_after == 0) break;
        
        //Orthocenter
        num_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles inserting ORTHOCENTER: " <<num_obtuses_before<<endl;
        /*Steiner Points Orthocenter*/
        insert_orthocenter(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles after inserting ORTHOCENTER: "<<num_obtuses_after<<endl;
        //CGAL::draw(custom_cdt);
        if(num_obtuses_after == 0) break;

        //Incenter
        num_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles inserting INCENTER: " <<num_obtuses_before<<endl;
        /*Steiner Points Incenter*/
        insert_incenter(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles after inserting INCENTER: "<<num_obtuses_after<<endl;
        //CGAL::draw(custom_cdt);
        if(num_obtuses_after == 0) break;

        //Bisector
        num_obtuses_before = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles inserting BISECTOR: " <<num_obtuses_before<<endl;
        /*Steiner Points Bisector*/
        insert_bisector(custom_cdt, polygon);
        num_obtuses_after = count_obtuse_triangles(custom_cdt, polygon);
        cout<<"Number of obtuse triangles after inserting BISECTOR: "<<num_obtuses_after<<endl;
        CGAL::draw(custom_cdt);
        
        int end = num_obtuses_after;
        if(end < start ) progress = true;
        else progress = false;
        if(end == 0) break;
    }
    return 0;
}
