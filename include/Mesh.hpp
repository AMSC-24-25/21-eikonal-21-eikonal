#ifndef MESH_HPP
#define MESH_HPP

#include <array>
#include <vector>
#include <memory>
#include <fstream>
#include <iostream>
#include <string>
#include <stdexcept>

#include "Eikonal_traits.hpp"
#include "Node.hpp"
#include "MeshElement.hpp"

template<std::size_t PHDIM>
class Mesh {
public:
    using Point = typename Eikonal::Eikonal_traits<PHDIM>::Point;
    using NodePtr = std::shared_ptr<Node<PHDIM>>;

    void resize(const std::size_t size) {
        mesh_elements.reserve(size);
    }

    std::vector<Point> Points;
    std::vector<NodePtr> nodes;
    std::vector<Mesh_element<PHDIM>> mesh_elements;
};

#endif