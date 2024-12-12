#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

#include "solveEikonalLocalProblem.hpp"

const double INF = 10e7;
const double EPSILON = 1e-6;
constexpr std::size_t PHDIM = 2;

struct Node
{
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

    unsigned int id;
    double u;
    bool isSource;
    Point p;
};

struct Mesh_element
{
    std::array<Node *, PHDIM + 1u> vertex;
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
    int k = 0;
    while (!activeList.empty() && k < 1000)
    {
        k++;
        std::vector<int> toAdd;
        std::vector<int> toRemove;

        for (auto it = activeList.begin(); it != activeList.end(); ++it)
        {
            Node *node = nodes[*it];
            double previous_value = node->u;
            node->u = solveLocal(*node);

            if (std::abs(previous_value - node->u) < EPSILON)
            {
                for (auto &neighbour : getNeighbours(*node))
                {
                    if (std::find(activeList.begin(), activeList.end(), neighbour->id) == activeList.end()
                        && !neighbour->isSource)
                    {
                        double p = neighbour->u;
                        double q = solveLocal(*neighbour);
                        if (p > q)
                        {
                            neighbour->u = q;
                            toAdd.push_back(neighbour->id);
                            std::cout << neighbour->id << std::endl;
                        }
                    }
                }
                toRemove.push_back(*it);
            }
        }

        for (const auto &id : toRemove)
        {
            activeList.erase(std::remove(activeList.begin(), activeList.end(), id), activeList.end());
        }

        for (const auto &id : toAdd)
        {
            activeList.push_back(id);
        }
    }
}

    void printResults() const
    {
        for (const auto &pair : nodes)
        {
            std::cout << "Node " << pair.second->id << ": u = " << pair.second->u << std::endl;
        }
    }

    std::vector<Node *> getNeighbours(Node &node)
    {
        std::unordered_set<unsigned int> neighbour_ids;
        std::vector<Node *> neighbours;

        for (auto &mesh_element : nodeToElements[node.id])
        {
            for (auto &mesh_node : mesh_element.vertex)
            {
                if (node.id != mesh_node->id && neighbour_ids.find(mesh_node->id) == neighbour_ids.end())
                {
                    neighbour_ids.insert(mesh_node->id);
                    neighbours.push_back(nodes[mesh_node->id]);
                }
            }
        }

        return neighbours;
    }

private:
    std::vector<Mesh_element> &mesh;
    std::unordered_map<unsigned int, Node *> nodes;
    std::unordered_map<unsigned int, std::vector<Mesh_element>> nodeToElements;
    std::vector<int> activeList;

    bool isInActiveList(Node &node)
    {
        return std::find(activeList.begin(), activeList.end(), node.id) != activeList.end();
    }

    void initializeMaps()
    {
        for (size_t i = 0; i < mesh.size(); ++i)
        {
            for (auto &node : mesh[i].vertex)
            {
                nodes[node->id] = node;
                nodeToElements[node->id].push_back(mesh[i]);
            }
        }
    }

    void initialize()
    {
        for (auto &m_element : mesh)
        {
            for (auto &node : m_element.vertex)
            {
                if (node->isSource)
                {
                    node->u = 0.0;
                    // activeList.push_back(node.id);
                    for(auto &neighbour : getNeighbours(*node)){
                        if(!isInActiveList(*neighbour) && !neighbour->isSource){
                            activeList.push_back(neighbour->id);
                        }
                    }
                }
                else
                {
                    node->u = INF;
                }
            }
        }
    }

    double solveLocal(Node &node)
    {
        using Point = Eikonal::Eikonal_traits<PHDIM>::Point;
        using VectorExt = Eikonal::Eikonal_traits<PHDIM>::VectorExt;

        double min_value = INF;
        std::vector<Node *> nodes_for_points;

        for (auto &mesh_element : nodeToElements[node.id])
        {
            for (auto &mesh_node : mesh_element.vertex)
            {
                if (mesh_node->id != node.id)
                {
                    nodes_for_points.push_back(nodes[mesh_node->id]);
                }
            }
            
            Eikonal::SimplexData<PHDIM> simplex{{nodes_for_points[0]->p, nodes_for_points[1]->p, node.p}};
            VectorExt values;
            values << nodes_for_points[0]->u, nodes_for_points[1]->u;

            Eikonal::solveEikonalLocalProblem<PHDIM> solver{simplex, values};
            auto sol = solver();
            
            min_value = std::min(min_value, sol.value);

            nodes_for_points.clear();
        }

        return min_value;
    }
};

int main()
{
    using Point = Eikonal::Eikonal_traits<PHDIM>::Point;

    Point p1, p2, p3, p4, p5, p6, p7, p8, p9;
    
    p1 << 0., 0.;
    p2 << 0., 1.;
    p3 << 1., 1.;
    p4 << 1., 0.;
    p5 << 0.0, 2.0;
    p6 << 2.0, 2.0;
    p7 << 2.0, 0.0;
    p8 << 1.0, 2.0;
    p9 << 2.0, 1.0;

    Node n1 = {0, 0.0, true, p1};
    Node n2 = {1, 0.0, true, p2};
    Node n3 = {2, 0.0, false, p3};
    Node n4 = {3, 0.0, true, p4};
    Node n5 = {4, 0.0, false, p5};
    Node n6 = {5, 0.0, false, p6};
    Node n7 = {6, 0.0, false, p7};
    Node n8 = {7, 0.0, false, p8};
    Node n9 = {8, 0.0, false, p9};

    Mesh_element m1 = {{&n1, &n2, &n3}};
    Mesh_element m2 = {{&n1, &n4, &n3}};
    Mesh_element m3 = {{&n2, &n3, &n5}};
    Mesh_element m4 = {{&n3, &n5, &n8}};
    Mesh_element m5 = {{&n4, &n7, &n3}};
    Mesh_element m6 = {{&n3, &n9, &n7}};
    Mesh_element m7 = {{&n3, &n9, &n8}};
    Mesh_element m8 = {{&n8, &n9, &n6}};

    std::vector<Mesh_element> mesh = {m1, m2, m3, m4, m5, m6, m7, m8};

    EikonalSolver solver(mesh);
    solver.printResults();
    solver.update();
    solver.printResults();

    return 0;
}

