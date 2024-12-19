#ifndef MESHELEMENT_HPP
#define MESHELEMENT_HPP

#include <array>
#include <memory>
#include "Node.hpp"

template<unsigned int PHDIM>
struct Mesh_element {
    std::array<NodePtr<PHDIM>, PHDIM + 1> vertex;
};

#endif // MESHELEMENT_HPP