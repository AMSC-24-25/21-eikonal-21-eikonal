#ifndef MESHELEMENT_HPP
#define MESHELEMENT_HPP

#include <array>
#include <memory>
#include "Node.hpp"

struct Mesh_element {
    std::array<NodePtr, 3> vertex;
};

#endif // MESHELEMENT_HPP