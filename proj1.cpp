#include "includes/utils/functions.h"

using namespace boost::json;
using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using CDT = CGAL::Constrained_Delaunay_triangulation_2<K>;
using Point = CDT::Point;
using Custom_CDT = Custom_Constrained_Delaunay_triangulation_2<K>;

int main() {
    value jv;
    read_json("data10.json", jv);

    std::vector<Point> points;
    std::vector<std::pair<int, int>> additional_constraints;
    std::vector<int> region_boundary;

    if (jv.is_object()) {
        const auto& obj = jv.as_object();
        const auto& x_array = obj.at("points_x").as_array();
        const auto& y_array = obj.at("points_y").as_array();
        const auto& boundary_array = obj.at("region_boundary").as_array();
        const auto& constraints_array = obj.at("additional_constraints").as_array();

        for (std::size_t i = 0; i < x_array.size(); ++i) {
            double x = x_array[i].is_double() ? x_array[i].as_double() : static_cast<double>(x_array[i].as_int64());
            double y = y_array[i].is_double() ? y_array[i].as_double() : static_cast<double>(y_array[i].as_int64());
            points.emplace_back(x, y);
        }

        for (const auto& idx : boundary_array) {
            region_boundary.push_back(idx.as_int64());
        }
        //std::vector<Point> region_boundary_edges = make_region_boundary(region_boundary, points);

        for (const auto& constraint : constraints_array) {
            int idx1 = constraint.as_array()[0].as_int64();
            int idx2 = constraint.as_array()[1].as_int64();
            if (idx1 < points.size() && idx2 < points.size()) {
                additional_constraints.emplace_back(idx1, idx2);
            }
        }
    }

    //CDT cdt;
    Custom_CDT custom_cdt;
    std::vector<CDT::Vertex_handle> vertex_handles;
    for (const auto& point : points) {
        vertex_handles.push_back(custom_cdt.insert(point));
    }

    
    for (size_t i = 0; i < region_boundary.size(); ++i) {
        int idx1 = region_boundary[i];
        int idx2 = region_boundary[(i + 1) % region_boundary.size()];
        custom_cdt.insert_constraint(points[idx1], points[idx2]);
    }

    for (const auto& constraint : additional_constraints) {
        custom_cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }
    
    //CGAL::draw(cdt);
    int num_obtuses_before = count_obtuse_triangles(custom_cdt);
    std::cout << "Number of obtuse triangles before flips: " << num_obtuses_before << std::endl;

    /*Flips*/
    start_the_flips(custom_cdt, points, additional_constraints, region_boundary);
    int num_obtuses_after = count_obtuse_triangles(custom_cdt);
    std::cout << "Number of obtuse triangles after flips cdt: " << num_obtuses_after << std::endl;
    CGAL::draw(custom_cdt);
    /*==========================================================================*/

    /*Steiner Points*/
    //Custom_CDT custom_cdt;
    /*for (auto vit = cdt.finite_vertices_begin(); vit != cdt.finite_vertices_end(); ++vit) {
        custom_cdt.insert_no_flip(vit->point());
    }
     // Check for validity
    if (!custom_cdt.is_valid()) {
        std::cerr << "Warning: Triangulation is not valid after copying." << std::endl;
    }

    int num_obtuses_before_steiner = count_obtuse_triangles(custom_cdt);
    std::cout << "Number of obtuse triangles before inserting Steiner points custom_cdt: " << num_obtuses_before_steiner << std::endl;
    CGAL::draw(custom_cdt);

    for (const auto& constraint : additional_constraints) {
        custom_cdt.insert_constraint(points[constraint.first], points[constraint.second]);
    }
    
    CGAL::draw(custom_cdt);*/

    insert_steiner_points(custom_cdt, points, region_boundary);

    int num_obtuses_after_steiner = count_obtuse_triangles(custom_cdt);
    std::cout << "Number of obtuse triangles after inserting Steiner points custom_cdt: " << num_obtuses_after_steiner << std::endl;

    CGAL::draw(custom_cdt);

    return 0;
}
