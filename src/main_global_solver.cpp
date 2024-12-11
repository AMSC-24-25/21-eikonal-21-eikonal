#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <omp.h>
#include <unordered_map>

#include "solveEikonalLocalProblem.hpp"

const double INF = std::numeric_limits<double>::infinity();
const double EPSILON = 1e-6;
constexpr std::size_t PHDIM = 2;

struct Node {
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

    unsigned int id;
    double u;
    bool isSource;
    Point p;
    //std::vector<int> neighbors;
};

struct Mesh_element {
    std::array<Node,PHDIM+1u> vertex;
};

class EikonalSolver {
public:
    EikonalSolver(std::vector<Mesh_element>& mesh) : mesh(mesh) {
        initialize();
    }

    void solve() {
        update();
    }

    void printResults() const {
        for (const auto& node : nodes) {
            std::cout << "Node " << node.second.id << ": u = " << node.second.u << std::endl;
        }
    }

private:
    std::vector<Mesh_element>& mesh;
    std::unordered_map<unsigned int, Node> nodes;
    std::unordered_map<unsigned int, std::vector<Mesh_element>> nodeToElements;
    std::queue<int> activeList;

    void initialize() {
        for (auto& m_element : mesh) {
            for (auto& node : m_element.vertex) {
                if (node.isSource) {
                    node.u = 0.0;
                } else {
                    node.u = INF;
                }
            }
        }

        for (const auto& pair : nodeToElements) {
            for (auto& mesh_element : pair.second) {
                for (auto& node : mesh_element.vertex) {
                    if (pair.first != node.id && node.isSource) {
                        activeList.push(pair.first);
                    }
                }
            }
        }
    }

    double solveLocalProblem(const Node& node) {
        // Implement the local solver here
        // This is a placeholder implementation
        return node.u;
    }

    void update() {
        while (!activeList.empty()) {
            int currentNodeId = activeList.front();
            activeList.pop();

            Node& currentNode = nodes[currentNodeId];
            double oldU = currentNode.u;
            double newU = solveLocalProblem(currentNode);

            if (std::abs(oldU - newU) < EPSILON) {
                currentNode.u = newU;
                for (int neighborId : currentNode.neighbors) {
                    Node& neighbor = nodes[neighborId];
                    if (neighbor.u > newU) {
                        neighbor.u = newU;
                        activeList.push(neighborId);
                    }
                }
            }
        }
    }
};

int main() {
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;
    
    // Example setup for 2D triangular mesh or 3D tetrahedral mesh
    
    Point p1,p2,p3,p4;
    p1<<0.,0.;
    p2<<0.,1.;
    p3<<1.,1.;
    p4<<1.,0.;

    Node n1, n2, n3, n4;
    n1 = {0, 0.0, true, p1};
    n2 = {1, 0.0, true, p2};
    n3 = {2, 0.0, false, p3};
    n4 = {3, 0.0, true, p4};

    Mesh_element m1, m2;
    m1.vertex = {n1, n2, n3};
    m2.vertex = {n1, n4, n3};

    std::vector<Mesh_element> mesh = {m1, m2};

    EikonalSolver solver(mesh);
    solver.solve();
    solver.printResults();

    return 0;
}
