#ifndef LOAD_MESH_HPP 
#define LOAD_MESH_HPP

#include <array>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>
#include "Mesh.hpp"

template<std::size_t PHDIM>
class loadMesh {
public:
    static std::vector<Mesh_element<PHDIM>> init_Mesh(const std::string& mesh_path, Mesh<PHDIM>& mesh) {
        std::ifstream mesh_file(mesh_path);
        if (!mesh_file.is_open()) {
            throw std::runtime_error("Unable to open mesh file: " + mesh_path);
        }

        mesh.Points.clear();
        mesh.nodes.clear();
        mesh.mesh_elements.clear();

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
        std::string data_type;
        mesh_file >> data_type; // This captures "float" or "double"

        mesh.Points.reserve(num_of_vertices);
        mesh.nodes.reserve(num_of_vertices);
        for (int i = 0; i < num_of_vertices; ++i) {
            typename Mesh<PHDIM>::Point p;
            for (std::size_t dim = 0; dim < PHDIM; ++dim) {
                if (!(mesh_file >> p[dim])) {
                    throw std::runtime_error("Error reading vertex coordinates");
                }
            }
            if (PHDIM == 2) {
                double ignore;
                mesh_file >> ignore; // Ignore the z-coordinate for 2D meshes
            }

            mesh.Points.push_back(p);
            auto node = std::make_shared<Node<PHDIM>>(Node<PHDIM>{
                static_cast<unsigned int>(i),
                0.0,
                false,
                p
            });
            mesh.nodes.push_back(node);
        }

        mesh_file >> section_marker;
        if (section_marker != "POLYGONS") {
            throw std::runtime_error("Expected POLYGONS section in VTK file");
        }

        int num_of_mesh_elements;
        int total_indices;
        mesh_file >> num_of_mesh_elements >> total_indices;

        mesh.mesh_elements.reserve(num_of_mesh_elements);
        for (int i = 0; i < num_of_mesh_elements; ++i) {
            int vertices_per_element;
            mesh_file >> vertices_per_element;

            if (vertices_per_element != PHDIM + 1) {
                throw std::runtime_error("Invalid number of vertices per element");
            }

            std::array<NodePtr<PHDIM>, PHDIM + 1> element_nodes;
            for (std::size_t j = 0; j < vertices_per_element; ++j) {
                unsigned int vertex_index;
                mesh_file >> vertex_index;

                if (vertex_index >= mesh.nodes.size()) {
                    throw std::runtime_error("Vertex index out of bounds");
                }

                element_nodes[j] = mesh.nodes[vertex_index];
            }

            mesh.mesh_elements.push_back(Mesh_element<PHDIM>{element_nodes});
        }

        mesh_file.close();
        return mesh.mesh_elements;
    }
};

// Explicit instantiation for common use cases
template class loadMesh<2>;
template class loadMesh<3>;

#endif // LOAD_MESH_HPP