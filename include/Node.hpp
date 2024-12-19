#ifndef NODE_HPP
#define NODE_HPP

#include "Eikonal_traits.hpp"
#include <memory>

template<unsigned int PHDIM>
struct Node
{
    using Point = typename Eikonal::Eikonal_traits<PHDIM>::Point;

    unsigned int id;
    double u;
    bool isSource;
    Point p;
};

template<unsigned int PHDIM>
using NodePtr = std::shared_ptr<Node<PHDIM>>;

#endif // NODE_HPP