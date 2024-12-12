#ifndef EIKONALSOLVER_HPP
#define EIKONALSOLVER_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "MeshElement.hpp"

class EikonalSolver
{
public:
    EikonalSolver(std::vector<Mesh_element> &mesh);
    void update();
    void printResults() const;
    std::vector<NodePtr> getNeighbours(Node &node);

private:
    std::vector<Mesh_element> &mesh;
    std::unordered_map<unsigned int, NodePtr> nodes;
    std::unordered_map<unsigned int, std::vector<Mesh_element>> nodeToElements;
    std::vector<int> activeList;

    bool isInActiveList(Node &node);
    void initializeMaps();
    void initialize();
    double solveLocal(Node &node);
};

#endif // EIKONALSOLVER_HPP