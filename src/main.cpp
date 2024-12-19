#include <iostream>
#include <vector>
#include <memory>
// Remove comments from the following line to enable parallel solver
//#include "ParallelEikonalSolver.hpp"
#include "EikonalSolver.hpp"
#include <chrono>
#include "Node.hpp"
#include "MeshElement.hpp"
#include "Mesh.hpp"
#include "loadMesh.hpp"
#include "VTKWriter.hpp"


int main()
{
#if DIMENSION == 2
    constexpr unsigned int PHDIM = 2;
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;
    using Mat = typename Eikonal::Eikonal_traits<PHDIM>::MMatrix;

    // Create mesh object
    Mesh<PHDIM> mesh;

    // Load mesh from file
    try
    {
        loadMesh<PHDIM>::init_Mesh("../tests/mesh2D.vtk", mesh);
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
    M_matrix << 1.0, 0.0,
        0.0, 1.0;

#else
    constexpr unsigned int PHDIM = 3;
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;
    using Mat = typename Eikonal::Eikonal_traits<PHDIM>::MMatrix;

    // Create mesh object
    Mesh<PHDIM> mesh;

    // Load mesh from file
    try
    {
        loadMesh<PHDIM>::init_Mesh("../tests/mesh3D.vtk", mesh);
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
    M_matrix << 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 1.0;

#endif

    // Initialize and run solver
    // ParallelEikonalSolver<PHDIM> solver(mesh.mesh_elements, M_matrix);
    EikonalSolver<PHDIM> solver(mesh.mesh_elements, M_matrix);
    solver.printResults();

    // Measure time for update
    // auto start = std::chrono::high_resolution_clock::now();
    // solver.update();
    // auto end = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> duration = end - start;
    // std::cout << "Time taken for update: " << duration.count() << " seconds" << std::endl;

    solver.update();

    solver.printResults();

    // Write solution to VTK file
    try
    {
        VTKWriter<PHDIM>::write("solution.vtk", mesh);
        //VTKWriter<PHDIM>::write("parallel_solution.vtk", mesh);
    }
    catch (const std::runtime_error &e)
    {
        std::cerr << "Error writing VTK file: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}