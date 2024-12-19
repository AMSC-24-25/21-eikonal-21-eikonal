#ifndef PARALLELEIKONALSOLVER_HPP
#define PARALLELEIKONALSOLVER_HPP

#include "EikonalSolver.hpp"
#include "MeshElement.hpp"
#include "solveEikonalLocalProblem.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <omp.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Parallel solver for the Eikonal equation using OpenMP.
 *
 * @tparam PHDIM Dimension of the problem space.
 */
template <unsigned int PHDIM> class ParallelEikonalSolver {
  using Mat = typename Eikonal::Eikonal_traits<PHDIM>::MMatrix;

public:
  /**
   * @brief Construct a new Parallel Eikonal Solver.
   *
   * @param mesh Reference to the mesh elements.
   * @param matrix Reference to the matrix used in calculations.
   */
  ParallelEikonalSolver(std::vector<Mesh_element<PHDIM>> &mesh, Mat &matrix)
      : mesh(mesh), mat(matrix) {
    initializeMaps();
    initialize();
  }

  /**
   * @brief Update the solution of the Eikonal equation.
   */
 void update()
    {
        while (!activeList.empty())
        {
            std::vector<int> toAdd;
            std::vector<int> toRemove;

            #pragma omp parallel default(shared)
            {
                std::vector<int> local_toAdd;
                std::vector<int> local_toRemove;

                #pragma omp for schedule(dynamic) nowait
                for (size_t idx = 0; idx < activeList.size(); ++idx)
                {
                    int node_id = activeList[idx];
                    NodePtr<PHDIM> node = nodes[node_id];

                    double previous_value = node->u;
                    double new_u = solveLocal(*node);
                    
                    #pragma omp atomic write
                    node->u = new_u;

                    if (std::abs(previous_value - new_u) < EPSILON)
                    {
                        auto neighbours = getNeighbours(*node);
                        for (auto& neighbour : neighbours)
                        {
                            if (!neighbour->isSource)
                            {
                                bool is_active = false;
                            #pragma omp critical
                            {
                                is_active = isInActiveList(*neighbour);
                            }
                            if (!is_active && !neighbour->isSource)
                            {
                                double p = neighbour->u;
                                double q = solveLocal(*neighbour);

                                    if (p > q)
                                    {
                                        #pragma omp atomic write
                                        neighbour->u = q;
                                        local_toAdd.push_back(neighbour->id);
                                    }
                                }
                            }
                        }
                        local_toRemove.push_back(node_id);
                    }
                }

                #pragma omp critical
                {
                    toAdd.insert(toAdd.end(), local_toAdd.begin(), local_toAdd.end());
                    toRemove.insert(toRemove.end(), local_toRemove.begin(), local_toRemove.end());
                }
            }

            #pragma omp single
            {
                for (const auto& id : toRemove)
                {
                    activeList.erase(std::remove(activeList.begin(), activeList.end(), id), activeList.end());
                }
                for (const auto& id : toAdd)
                {
                    activeList.push_back(id);
                }
            }
        }
    }

  /**
   * @brief Print the results of the computation.
   */
  void printResults() const {
    for (const auto &pair : nodes) {
      std::cout << "Node " << pair.second->id << ": u = " << pair.second->u
                << std::endl;
    }
  }

  /**
   * @brief Get the neighbours of a given node.
   *
   * @param node The node whose neighbours are to be found.
   * @return std::vector<NodePtr<PHDIM>> Vector of neighbouring node pointers.
   */
  std::vector<NodePtr<PHDIM>> getNeighbours(Node<PHDIM> &node) {
    std::unordered_set<unsigned int> neighbour_ids;
    std::vector<NodePtr<PHDIM>> neighbours;
    for (auto &mesh_element : nodeToElements[node.id]) {
      for (auto &mesh_node : mesh_element.vertex) {
        if (node.id != mesh_node->id &&
            neighbour_ids.find(mesh_node->id) == neighbour_ids.end()) {
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
  std::unordered_map<unsigned int, std::vector<Mesh_element<PHDIM>>>
      nodeToElements;
  std::vector<int> activeList;
  Mat &mat;

  bool isInActiveList(Node<PHDIM> &node) {
    return std::find(activeList.begin(), activeList.end(), node.id) !=
           activeList.end();
  }

  void initializeMaps() {
#pragma omp parallel for schedule(dynamic) default(shared)
    for (size_t i = 0; i < mesh.size(); ++i) {
#pragma omp critical
      for (auto &node : mesh[i].vertex) {
        nodes[node->id] = node;
        nodeToElements[node->id].push_back(mesh[i]);
      }
    }
  }

  void initialize() {
#pragma omp parallel for schedule(dynamic) default(shared)
    for (auto &m_element : mesh) {
      for (auto &node : m_element.vertex) {
        if (node->isSource) {
          node->u = 0.0;
#pragma omp critical
          for (auto &neighbour : getNeighbours(*node)) {
            if (!isInActiveList(*neighbour) && !neighbour->isSource) {
              activeList.push_back(neighbour->id);
            }
          }
        } else {
          node->u = INF;
        }
      }
    }
  }

  double solveLocal(Node<PHDIM> &node) {
    using Point = typename Eikonal::Eikonal_traits<PHDIM>::Point;
    using VectorExt = typename Eikonal::Eikonal_traits<PHDIM>::VectorExt;

    double min_value = INF;
    auto &elements = nodeToElements[node.id];

#pragma omp parallel default(shared) reduction(min : min_value)
    {
      std::vector<NodePtr<PHDIM>> nodes_for_points;

#pragma omp for schedule(dynamic) nowait
      for (size_t idx = 0; idx < elements.size(); ++idx) {
        auto &mesh_element = elements[idx];

        for (auto &mesh_node : mesh_element.vertex) {
          if (mesh_node->id != node.id) {
            nodes_for_points.push_back(mesh_node);
          }
        }

        if (nodes_for_points.size() < PHDIM) {
          nodes_for_points.clear();
          continue;
        }

        std::array<Point, PHDIM + 1> simplex_points;
        for (unsigned int i = 0; i < PHDIM; ++i) {
          simplex_points[i] = nodes_for_points[i]->p;
        }
        simplex_points[PHDIM] = node.p;

        VectorExt values;
        for (unsigned int i = 0; i < PHDIM; ++i) {
          values[i] = nodes_for_points[i]->u;
        }

        Eikonal::SimplexData<PHDIM> simplex{simplex_points, mat};
        Eikonal::solveEikonalLocalProblem<PHDIM> solver{simplex, values};
        auto sol = solver();
        min_value = std::min(min_value, sol.value);
        nodes_for_points.clear();
      }
    }
    return min_value;
  }
};

#endif // PARALLELEIKONALSOLVER_HPP