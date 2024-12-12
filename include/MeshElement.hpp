
#ifndef MESHELEMENT_HPP
#define MESHELEMENT_HPP

#include <array>
#include "Node.hpp"

struct Mesh_element
{
    std::array<Node *, 3> vertex;
};

#endif // MESHELEMENT_HPP