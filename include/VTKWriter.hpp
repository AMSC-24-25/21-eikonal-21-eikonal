#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include <fstream>
#include <string>
#include "Mesh.hpp"

template<unsigned int PHDIM>
class VTKWriter {
public:
    static void write(const std::string& filename, const Mesh<PHDIM>& mesh) {
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Unable to open file for writing: " + filename);
        }

        // Write VTK header
        file << "# vtk DataFile Version 3.0\n";
        file << "Eikonal solution\n";
        file << "ASCII\n";
        file << "DATASET UNSTRUCTURED_GRID\n";

        // Write points
        file << "POINTS " << mesh.nodes.size() << " double\n";
        for (const auto& node : mesh.nodes) {
            for (size_t i = 0; i < PHDIM; ++i) {
                file << node->p[i] << " ";
            }
            if (PHDIM == 2) file << "0 "; // Add z=0 for 2D
            file << "\n";
        }

        // Write cells
        const size_t cells_size = mesh.mesh_elements.size();
        const size_t points_per_cell = PHDIM + 1;
        file << "\nCELLS " << cells_size << " " << cells_size * (points_per_cell + 1) << "\n";
        for (const auto& element : mesh.mesh_elements) {
            file << points_per_cell << " ";
            for (size_t i = 0; i < points_per_cell; ++i) {
                file << element.vertex[i]->id << " ";
            }
            file << "\n";
        }

        // Write cell types
        file << "\nCELL_TYPES " << cells_size << "\n";
        for (size_t i = 0; i < cells_size; ++i) {
            file << (PHDIM == 2 ? 5 : 10) << "\n"; // 5 for triangle, 10 for tetrahedron
        }

        // Write point data (solution values)
        file << "\nPOINT_DATA " << mesh.nodes.size() << "\n";
        file << "SCALARS solution double 1\n";
        file << "LOOKUP_TABLE default\n";
        for (const auto& node : mesh.nodes) {
            file << node->u << "\n";
        }

        file.close();
    }
};

#endif // VTK_WRITER_HPP
