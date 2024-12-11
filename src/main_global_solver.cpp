#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>

#include "solveEikonalLocalProblem.hpp"

const double INF = std::numeric_limits<double>::infinity();
const double EPSILON = 1e-6;
constexpr std::size_t PHDIM = 2;

struct Node
{
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

    unsigned int id;
    double u;
    bool isSource;
    Point p;
    // std::vector<int> neighbors;
};

struct Mesh_element
{
    std::array<Node, PHDIM + 1u> vertex;
};

class EikonalSolver
{
public:
    EikonalSolver(std::vector<Mesh_element> &mesh) : mesh(mesh)
    {
        initializeMaps();
        initialize();
    }

    void update()
    {
        int k=0;
        while (!activeList.empty() && k<10)
        {
            k++;
            std::cout << activeList.size() << std::endl;
            for (auto &nodeId : activeList)
            {
                auto &node = nodes[nodeId];
                double previous_value = node.u;
                node.u = solveLocal(node);

                if (abs(previous_value - node.u) < EPSILON)
                {
                    for (auto &neighbourg : getNeighbourgs(node))
                    {
                        if (std::find(activeList.begin(), activeList.end(), neighbourg.id) == activeList.end()) {
                            auto &p = neighbourg.u;
                            auto q = solveLocal(neighbourg);
                            if (p > q) {
                                neighbourg.u = q;
                                activeList.push_back(neighbourg.id);
                            }
                        }
                    }
                    activeList.erase(std::remove(activeList.begin(), activeList.end(), nodeId), activeList.end());
                }
            }
        }
    }

    void printResults() const
    {
        for (const auto &node : nodes)
        {
            std::cout << "Node " << node.second.id << ": u = " << node.second.u << std::endl;
        }
    }

    std::vector<Node> getNeighbourgs(Node &node)
    {
        std::unordered_set<unsigned int> neighbourg_ids;
        std::vector<Node> neighbourgs;

        for (auto &mesh_element : nodeToElements[node.id])
        {
            for (auto &mesh_node : mesh_element.vertex)
            {
                if (node.id != mesh_node.id && neighbourg_ids.find(mesh_node.id) == neighbourg_ids.end())
                {
                    neighbourg_ids.insert(mesh_node.id);
                    neighbourgs.push_back(mesh_node);
                }
            }
        }

        return neighbourgs;
    }

//private:
    std::vector<Mesh_element> &mesh;
    std::unordered_map<unsigned int, Node> nodes;
    std::unordered_map<unsigned int, std::vector<Mesh_element>> nodeToElements;
    std::vector<int> activeList;

    void initializeMaps()
    {
        for (size_t i = 0; i < mesh.size(); ++i)
        {
            for (const auto &node : mesh[i].vertex)
            {
                nodes.insert({node.id, node});
                nodeToElements[node.id].push_back(mesh[i]);
            }
        }
    }

    void initialize()
    {
        for (auto &m_element : mesh)
        {
            for (auto &node : m_element.vertex)
            {
                if (node.isSource)
                {
                    node.u = 0.0;
                }
                else
                {
                    node.u = INF;
                }
            }
        }

        for (const auto &pair : nodeToElements)
        {
            for (auto &neighbourg : getNeighbourgs(nodes[pair.first]))
            {
                if (nodes[pair.first].isSource)
                {
                    activeList.push_back(pair.first);
                }
            }
        }
    }

    double solveLocal(Node &node)
    {
        using Point = Eikonal::Eikonal_traits<PHDIM>::Point;
        using VectorExt = Eikonal::Eikonal_traits<PHDIM>::VectorExt;

        double min_value = INF;
        std::vector<Node> nodes_for_points;
        for (auto &mesh_element : nodeToElements[node.id])
        {
            for (auto &mesh_node : mesh_element.vertex)
            {
                if (mesh_node.id != node.id)
                {
                    nodes_for_points.push_back(mesh_node);
                }
            }

            Eikonal::SimplexData<PHDIM> simplex{{nodes_for_points[0].p, nodes_for_points[1].p, node.p}};
            VectorExt values;
            values << nodes_for_points[0].u, nodes_for_points[1].u;

            Eikonal::solveEikonalLocalProblem<PHDIM> solver{simplex, values};
            auto sol = solver();
            double min_value = std::min(min_value, sol.value);
            
            nodes_for_points.clear();
        }
        return min_value;
    }
};

int main()
{
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

    // Example setup for 2D triangular mesh or 3D tetrahedral mesh

    Point p1, p2, p3, p4;
    p1 << 0., 0.;
    p2 << 0., 1.;
    p3 << 1., 1.;
    p4 << 1., 0.;

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
    solver.printResults();
    solver.update();
    solver.printResults();

    
    return 0;
}
