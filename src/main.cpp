#include <iostream>
#include <vector>
#include <memory>
#include "EikonalSolver.hpp"
#include "Node.hpp"
#include "MeshElement.hpp"
#include "Mesh.hpp"
#include "loadMesh.hpp"

int main()
{
    constexpr unsigned int PHDIM = 2;
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

    // Create mesh object
    Mesh<PHDIM> mesh;
    
    // Load mesh from file
    try {
        loadMesh<PHDIM>::init_Mesh("../tests/m.vtk", mesh);
    } catch (const std::runtime_error& e) {
        std::cerr << "Error loading mesh: " << e.what() << std::endl;
        return 1;
    }

    // Set source nodes (example: set first node as source)
    if (!mesh.nodes.empty()) {
        mesh.nodes[1]->isSource = true;
    }

    // Initialize and run solver
    EikonalSolver<PHDIM> solver(mesh.mesh_elements);
    solver.printResults();
    solver.update();
    solver.printResults();

    return 0;
}