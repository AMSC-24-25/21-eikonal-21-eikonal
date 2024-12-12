#ifndef NODE_HPP
#define NODE_HPP

#include "Eikonal_traits.hpp"
#include <memory>

struct Node
{
    using Point = Eikonal::Eikonal_traits<2>::Point;

    unsigned int id;
    double u;
    bool isSource;
    Point p;
};

using NodePtr = std::shared_ptr<Node>;

#endif // NODE_HPP