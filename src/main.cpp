#include <iostream>
#include <vector>
#include <memory>
#include "EikonalSolver.hpp"
#include "Node.hpp"
#include "MeshElement.hpp"
#include "Mesh.hpp"
#include "loadMesh.hpp"
#include "VTKWriter.hpp"

int main()
{
#if DIMENSION == 2
    constexpr unsigned int PHDIM = 2;
#else
    constexpr unsigned int PHDIM = 3;
#endif

    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;
    using Mat = typename Eikonal::Eikonal_traits<PHDIM>::MMatrix;

    // Create mesh object
    Mesh<PHDIM> mesh;

    // Load mesh from file
    try
    {
#if DIMENSION == 2
        loadMesh<PHDIM>::init_Mesh("../tests/mesh2D.vtk", mesh);
#else
        loadMesh<PHDIM>::init_Mesh("../tests/mesh3D.vtk", mesh);
#endif
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error loading mesh: " << e.what() << std::endl;
        return 1;
    }

    // Set source nodes (example: set first node as source)
    if (!mesh.nodes.empty())
    {
        mesh.nodes[51 - 7]->isSource = true;
    }

    // Create anisotropy matrix
    Mat M_matrix;
#if DIMENSION == 2
    M_matrix << 1.0, 0.0,
        0.0, 1.0;
#else
    M_matrix << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;
#endif

    // Initialize and run solver
    EikonalSolver<PHDIM> solver(mesh.mesh_elements, M_matrix);
    solver.printResults();
    solver.update();
    solver.printResults();

    // Write solution to VTK file
    try
    {
        VTKWriter<PHDIM>::write("solution.vtk", mesh);
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error writing VTK file: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}