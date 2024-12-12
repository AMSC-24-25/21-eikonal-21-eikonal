#ifndef EIKONALSOLVER_HPP
#define EIKONALSOLVER_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>
#include "MeshElement.hpp"
#include "EikonalSolver.hpp"
#include <iostream>
#include <algorithm>
#include <cmath>
#include "solveEikonalLocalProblem.hpp"

const double INF = 10e7;
const double EPSILON = 1e-6;

template<unsigned int PHDIM>
class EikonalSolver
{
using Mat = typename Eikonal::Eikonal_traits<PHDIM>::MMatrix;
public:
    EikonalSolver(std::vector<Mesh_element<PHDIM>> &mesh, Mat &matrix) : mesh(mesh), mat(matrix)
    {
        initializeMaps();
        initialize();
    }

    void update()
    {
        while (!activeList.empty())
        {
            std::vector<int> toAdd;
            std::vector<int> toRemove;

            for (auto it = activeList.begin(); it != activeList.end(); ++it)
            {
                NodePtr<PHDIM> node = nodes[*it];
                double previous_value = node->u;
                node->u = solveLocal(*node);

                if (std::abs(previous_value - node->u) < EPSILON)
                {
                    for (auto &neighbour : getNeighbours(*node))
                    {
                        if (std::find(activeList.begin(), activeList.end(), neighbour->id) == activeList.end() && !neighbour->isSource)
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
    std::vector<NodePtr<PHDIM>> getNeighbours(Node<PHDIM> &node)
    {
        std::unordered_set<unsigned int> neighbour_ids;
        std::vector<NodePtr<PHDIM>> neighbours;

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
    std::vector<Mesh_element<PHDIM>> &mesh;
    std::unordered_map<unsigned int, NodePtr<PHDIM>> nodes;
    std::unordered_map<unsigned int, std::vector<Mesh_element<PHDIM>>> nodeToElements;
    std::vector<int> activeList;
    Mat& mat;

    bool isInActiveList(Node<PHDIM> &node)
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

    double solveLocal(Node<PHDIM> &node)
    {
        using Point = typename Eikonal::Eikonal_traits<PHDIM>::Point;
        using VectorExt = typename Eikonal::Eikonal_traits<PHDIM>::VectorExt;

        double min_value = INF;
        std::vector<NodePtr<PHDIM>> nodes_for_points;

        for (auto &mesh_element : nodeToElements[node.id])
        {
            for (auto &mesh_node : mesh_element.vertex)
            {
                if (mesh_node->id != node.id)
                {
                    nodes_for_points.push_back(nodes[mesh_node->id]);
                }
            }
            
            Eikonal::SimplexData<PHDIM> simplex{{nodes_for_points[0]->p, nodes_for_points[1]->p, node.p}, mat};
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

#endif // EIKONALSOLVER_HPP