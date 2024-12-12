
#ifndef NODE_HPP
#define NODE_HPP

#include "Eikonal_traits.hpp"

struct Node
{
    using Point = Eikonal::Eikonal_traits<2>::Point;

    unsigned int id;
    double u;
    bool isSource;
    Point p;
};

#endif // NODE_HPP