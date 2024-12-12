
#include "EikonalSolver.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include "solveEikonalLocalProblem.hpp"

const double INF = 10e7;
const double EPSILON = 1e-6;

EikonalSolver::EikonalSolver(std::vector<Mesh_element> &mesh) : mesh(mesh)
{
    initializeMaps();
    initialize();
}

void EikonalSolver::update()
{
    while (!activeList.empty())
    {
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

void EikonalSolver::printResults() const
{
    for (const auto &pair : nodes)
    {
        std::cout << "Node " << pair.second->id << ": u = " << pair.second->u << std::endl;
    }
}

std::vector<Node *> EikonalSolver::getNeighbours(Node &node)
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

bool EikonalSolver::isInActiveList(Node &node)
{
    return std::find(activeList.begin(), activeList.end(), node.id) != activeList.end();
}

void EikonalSolver::initializeMaps()
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

void EikonalSolver::initialize()
{
    for (auto &m_element : mesh)
    {
        for (auto &node : m_element.vertex)
        {
            if (node->isSource)
            {
                node->u = 0.0;
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

double EikonalSolver::solveLocal(Node &node)
{
    using Point = Eikonal::Eikonal_traits<2>::Point;
    using VectorExt = Eikonal::Eikonal_traits<2>::VectorExt;

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
        
        Eikonal::SimplexData<2> simplex{{nodes_for_points[0]->p, nodes_for_points[1]->p, node.p}};
        VectorExt values;
        values << nodes_for_points[0]->u, nodes_for_points[1]->u;

        Eikonal::solveEikonalLocalProblem<2> solver{simplex, values};
        auto sol = solver();
        
        min_value = std::min(min_value, sol.value);

        nodes_for_points.clear();
    }

    return min_value;
}