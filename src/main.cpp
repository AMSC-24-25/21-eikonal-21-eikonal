#include <iostream>
#include <vector>
#include <memory>
#include "EikonalSolver.hpp"
#include "Node.hpp"
#include "MeshElement.hpp"

int main()
{
    constexpr unsigned int PHDIM = 2;
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

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

    auto n1 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{0, 0.0, true, p1});
    auto n2 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{1, 0.0, true, p2});
    auto n3 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{2, 0.0, false, p3});
    auto n4 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{3, 0.0, true, p4});
    auto n5 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{4, 0.0, false, p5});
    auto n6 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{5, 0.0, false, p6});
    auto n7 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{6, 0.0, false, p7});
    auto n8 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{7, 0.0, false, p8});
    auto n9 = std::make_shared<Node<PHDIM>>(Node<PHDIM>{8, 0.0, false, p9});

    Mesh_element<PHDIM> m1 = {{n1, n2, n3}};
    Mesh_element<PHDIM> m2 = {{n1, n4, n3}};
    Mesh_element<PHDIM> m3 = {{n2, n3, n5}};
    Mesh_element<PHDIM> m4 = {{n3, n5, n8}};
    Mesh_element<PHDIM> m5 = {{n4, n7, n3}};
    Mesh_element<PHDIM> m6 = {{n3, n9, n7}};
    Mesh_element<PHDIM> m7 = {{n3, n9, n8}};
    Mesh_element<PHDIM> m8 = {{n8, n9, n6}};

    std::vector<Mesh_element<PHDIM>> mesh = {m1, m2, m3, m4, m5, m6, m7, m8};

    EikonalSolver<PHDIM> solver(mesh);
    solver.printResults();
    solver.update();
    solver.printResults();

    return 0;
}