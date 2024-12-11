#ifndef MESH_HPP
#define MESH_HPP

#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include "Eikonal_traits.hpp"

template<std::size_t PHDIM>
class Mesh {
public:
    using Point = typename Eikonal::Eikonal_traits<PHDIM>::Point;

    struct Node {
        unsigned int id;
        double u = 0.0;
        bool isSource = false;
        Point p;
    };

    struct Mesh_element {
        std::array<Node, PHDIM + 1> vertex;
    };

    Mesh() = default;

    std::vector<Mesh_element> init_Mesh(const std::string& mesh_path) {
        std::ifstream mesh_file(mesh_path);
        if (!mesh_file.is_open()) {
            throw std::runtime_error("Unable to open mesh file: " + mesh_path);
        }

        Points.clear();
        mesh_elements.clear();

        std::string header_line;
        for (int i = 0; i < 4; ++i) {
            if (!std::getline(mesh_file, header_line)) {
                throw std::runtime_error("Incomplete VTK header in mesh file");
            }
        }

        std::string section_marker;
        int num_of_vertices;
        mesh_file >> section_marker;
        if (section_marker != "POINTS") {
            throw std::runtime_error("Expected POINTS section in VTK file");
        }
        mesh_file >> num_of_vertices;
        mesh_file >> section_marker;

        Points.reserve(num_of_vertices);
        for (int i = 0; i < num_of_vertices; ++i) {
            Point p;
            for (std::size_t dim = 0; dim < PHDIM; ++dim) {
                if (!(mesh_file >> p[dim])) {
                    throw std::runtime_error("Error reading vertex coordinates");
                }
            }
            Points.push_back(p);
        }

        mesh_file >> section_marker;
        if (section_marker != "CELLS") {
            throw std::runtime_error("Expected CELLS section in VTK file");
        }
        
        int num_of_mesh_elements;
        mesh_file >> num_of_mesh_elements;
        mesh_file >> section_marker;

        mesh_elements.reserve(num_of_mesh_elements);
        for (int i = 0; i < num_of_mesh_elements; ++i) {
            int vertices_per_element;
            mesh_file >> vertices_per_element;
            
            /**if (vertices_per_element != static_cast<int>(PHDIM + 1)) {
                throw std::runtime_error("Vertex count mismatch for mesh dimension");
            }**/
            Mesh_element element;
            for (std::size_t j = 0; j < PHDIM + 1; ++j) {
                unsigned int vertex_index;
                mesh_file >> vertex_index;

                if (vertex_index >= Points.size()) {
                    throw std::runtime_error("Vertex index out of bounds");
                }

                element.vertex[j].id = vertex_index;
                element.vertex[j].p = Points[vertex_index];
            }
            mesh_elements.push_back(element);
        }

        mesh_file.close();
        return mesh_elements;
    }

    const std::vector<Point>& getPoints() const {
        return Points;
    }

    const std::vector<Mesh_element>& getMeshElements() const {
        return mesh_elements;
    }

private:
    const std::size_t D = PHDIM;
    std::vector<Point> Points;
    std::vector<Mesh_element> mesh_elements;
};

#endif