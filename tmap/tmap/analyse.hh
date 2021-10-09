/**
 * @file analyse.hh
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Analysing TMAP graphs and trees.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#ifndef ANALYSE_H
#define ANALYSE_H

#include <stdint.h>
#include <tuple>
#include <vector>
#include <stack>
#include <queue>
#include <utility>
#include <iostream>
#include <numeric>
#include <algorithm>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>

#include "layout.hh"

namespace tmap {
/**
 * @brief Get all clusters (uninterrupted).
 *
 * @param s Edge vertex a.
 * @param t Edge vertex b.
 * @param classes The classes according to which clusters are defined.
 * @param vectex_count The number of vertices in the graph.
 * @return std::vector<std::tuple<uint32_t, std::vector<uint32_t>>>
 */
std::vector<std::tuple<uint32_t, std::vector<uint32_t>>>
GetClusters(
  const tmap::GraphProperties& gp,
  const std::vector<uint32_t>& classes);

/**
 * @brief Hierarchical MST clustering (HEMST).
 *
 * @param gp A TMAP GraphProperties object.
 * @param k THe number of clusters.
 * @return std::vector<std::vector<uint32_t>>
 */
std::vector<std::vector<uint32_t>>
MSDR(tmap::GraphProperties gp);

}; // namespace tmap
#endif