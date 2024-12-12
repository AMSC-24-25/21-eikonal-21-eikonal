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
        std::array<Node, PHDIM+1 > vertex;
    };
    
    void resize(const std::size_t size) {
       
    }
    std::vector<Point> Points;

    std::vector<Mesh_element> mesh_elements;
};

#endif