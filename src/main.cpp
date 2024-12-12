
#include <iostream>
#include <vector>
#include "EikonalSolver.hpp"
#include "Node.hpp"
#include "MeshElement.hpp"

int main()
{
    using Point = Eikonal::Eikonal_traits<2>::Point;

    Point p1, p2, p3, p4, p5, p6, p7, p8, p9;
    
    p1 << 0., 0.;
    p2 << 0., 1.;
    p3 << 1., 1.;
    p4 << 1., 0.;
    p5 << 0.0, 2.0;
    p6 << 2.0, 2.0;
    p7 << 2.0, 0.0;
    p8 << 1.0, 2.0;
    p9 << 2.0, 1.0;

    Node n1 = {0, 0.0, true, p1};
    Node n2 = {1, 0.0, true, p2};
    Node n3 = {2, 0.0, false, p3};
    Node n4 = {3, 0.0, true, p4};
    Node n5 = {4, 0.0, false, p5};
    Node n6 = {5, 0.0, false, p6};
    Node n7 = {6, 0.0, false, p7};
    Node n8 = {7, 0.0, false, p8};
    Node n9 = {8, 0.0, false, p9};

    Mesh_element m1 = {{&n1, &n2, &n3}};
    Mesh_element m2 = {{&n1, &n4, &n3}};
    Mesh_element m3 = {{&n2, &n3, &n5}};
    Mesh_element m4 = {{&n3, &n5, &n8}};
    Mesh_element m5 = {{&n4, &n7, &n3}};
    Mesh_element m6 = {{&n3, &n9, &n7}};
    Mesh_element m7 = {{&n3, &n9, &n8}};
    Mesh_element m8 = {{&n8, &n9, &n6}};

    std::vector<Mesh_element> mesh = {m1, m2, m3, m4, m5, m6, m7, m8};

    EikonalSolver solver(mesh);
    solver.printResults();
    solver.update();
    solver.printResults();

    return 0;
}