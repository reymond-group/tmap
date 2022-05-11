/**
 * @file layout.hh
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Functions used for generating graph layouts from LSHForest instances
 * and edge lists.
 * @version 0.1
 * @date 2019-06-17
 *
 */

/* Comments for the documentation are (partially) taken from
   "An Experimental Evaluation of Multilevel Layout Methods"
   by G. Bartel et al. and the OGDF documentation. */

#ifndef LAYOUT_H
#define LAYOUT_H

#include <stdint.h>
#include <utility>
#include <tuple>
#include <vector>
#include <queue>
#include <algorithm>

#include <ogdf/basic/Graph.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/extended_graph_alg.h>
#include <ogdf/fileformats/GraphIO.h>

// #include <ogdf/energybased/FastMultipoleEmbedder.h>
#include <ogdf/energybased/multilevel_mixer/ScalingLayout.h>
#include <ogdf/energybased/multilevel_mixer/BarycenterPlacer.h>
#include <ogdf/energybased/multilevel_mixer/CirclePlacer.h>
#include <ogdf/energybased/multilevel_mixer/MedianPlacer.h>
#include <ogdf/energybased/multilevel_mixer/RandomPlacer.h>
#include <ogdf/energybased/multilevel_mixer/SolarPlacer.h>
#include <ogdf/energybased/multilevel_mixer/ZeroPlacer.h>

#include <ogdf/energybased/multilevel_mixer/EdgeCoverMerger.h>
#include <ogdf/energybased/multilevel_mixer/IndependentSetMerger.h>
#include <ogdf/energybased/multilevel_mixer/LocalBiconnectedMerger.h>
#include <ogdf/energybased/multilevel_mixer/SolarMerger.h>

#include <ogdf/basic/PreprocessorLayout.h>
#include <ogdf/packing/ComponentSplitterLayout.h>
#include <ogdf/packing/TileToRowsCCPacker.h>

#include <ogdf/energybased/DavidsonHarelLayout.h>
// #include <ogdf/energyb¦ased/DTreeMultilevelEmbedder.h>
#include <ogdf/energybased/FastMultipoleEmbedder.h>
#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/GEMLayout.h>
#include <ogdf/energybased/MultilevelLayout.h>
#include <ogdf/energybased/PivotMDS.h>
#include <ogdf/energybased/SpringEmbedderFRExact.h>
#include <ogdf/energybased/SpringEmbedderGridVariant.h>
#include <ogdf/energybased/SpringEmbedderKK.h>
#include <ogdf/energybased/StressMinimization.h>
#include <ogdf/energybased/TutteLayout.h>

#include <ogdf/energybased/NodeRespecterLayout.h>

#include <fstream>
#include <iostream>
#include <typeinfo>

namespace tmap
{
  // Forward declaration of LSHForest.
  class LSHForest;

  /**
   * @brief The placers available in OGDF.
   *
   */
  enum struct Placer
  {
    Barycenter =
        0,      ///< Places a vertex at the barycenter of its neighbors' position.
    Solar = 1,  ///< Uses information of the merging phase of the solar merger.
                ///< Places a new vertex on the direct line between two suns.
    Circle = 2, ///< Places the vertices in a circle around the barycenter and
                ///< outside of the current drawing
    Median = 3, ///< Places a vertex at the median position of the neighbor nodes
                ///< for each coordinate axis.
    Random = 4, ///< Places a vertex at a random position within the smallest
                ///< circle containing all vertices around the barycenter of the
                ///< current drawing.
    Zero = 5    ///< Places a vertex at the same position as its representative in
                ///< the previous level.
  };
  static const std::string placer_values[] = {"Barycenter", "Solar", "Circle",
                                              "Median", "Random", "Zero"};

  /**
   * @brief The mergers available in OGDF.
   *
   */
  enum struct Merger
  {
    EdgeCover =
        0, ///< Based on the matching merger. Computes an edge cover such that
           ///< each contained edge is incident to at least one unmatched vertex.
           ///< The cover edges are then used to merge their adjacent vertices.
    LocalBiconnected =
        1,     ///< Based on the edge cover merger. Avoids distortions by checking
               ///< whether biconnectivity will be lost in the local neighborhood
               ///< around the potential merging position.
    Solar = 2, ///< Vertices are partitioned into solar systems, consisting of
               ///< sun, planets and moons. The systems are then merged into the
               ///< sun vertices.
    IndependentSet =
        3 ///< Uses a maximal independent set filtration. See GRIP for details.
  };
  static const std::string merger_values[] = {"EdgeCover",
                                              "LocalBiconnected",
                                              "Solar",
                                              "IndependentSet"};

  /**
   * @brief The scaling types available in OGDF.
   *
   */
  enum struct ScalingType
  {
    Absolute = 0, ///< Absolute factor, can be used to scale relative to level
                  ///< size change.
    RelativeToAvgLength =
        1, ///< Scales by a factor relative to the average edge weights.
    RelativeToDesiredLength =
        2,                ///< Scales by a factor relative to the disired edge length.
    RelativeToDrawing = 3 ///< Scales by a factor relative to the drawing.
  };
  static const std::string scaling_types_values[] = {"Absolute",
                                                     "RelativeToAvgLength",
                                                     "RelativeToDesiredLength",
                                                     "RelativeToDrawing"};

  /**
   * @brief A struct containing all the configuration options available for and
   * applied to a layout.
   *
   */
  struct LayoutConfiguration
  {
    /**
     * @brief Construct a new Layout Configuration object.
     *
     */
    LayoutConfiguration()
        : k(10), kc(10), fme_iterations(1000), fme_randomize(false), fme_threads(4), fme_precision(4), sl_repeats(1), sl_extra_scaling_steps(2), sl_scaling_min(1.0), sl_scaling_max(1.0), sl_scaling_type(ScalingType::RelativeToDrawing), mmm_repeats(1), placer(Placer::Barycenter), merger(Merger::LocalBiconnected), merger_factor(2.0), merger_adjustment(0), node_size(1.0 / 65.0)
    {
    }

    /**
     * @brief Returns a string describing the set options.
     *
     * @return std::string
     */
    std::string ToString() const
    {
      return std::string("k: ") + std::to_string(k) + '\n' +
             "kc: " + std::to_string(kc) + '\n' +
             "fme_iterations: " + std::to_string(fme_iterations) + '\n' +
             "fme_randomize: " + std::to_string(fme_randomize) + '\n' +
             "fme_threads: " + std::to_string(fme_threads) + '\n' +
             "fme_precision: " + std::to_string(fme_threads) + '\n' +
             "sl_repeats: " + std::to_string(sl_repeats) + '\n' +
             "sl_extra_scaling_steps: " + std::to_string(sl_extra_scaling_steps) +
             '\n' + "sl_scaling_x: " + std::to_string(sl_scaling_min) + '\n' +
             "sl_scaling_y: " + std::to_string(sl_scaling_max) + '\n' +
             "sl_scaling_type: " + scaling_types_values[(int)sl_scaling_type] +
             '\n' + "mmm_repeats: " + std::to_string(mmm_repeats) + '\n' +
             "placer: " + placer_values[(int)placer] + '\n' +
             "merger: " + merger_values[(int)merger] + '\n' +
             "merger_factor: " + std::to_string(merger_factor) + '\n' +
             "merger_adjustment: " + std::to_string(merger_adjustment) + '\n' +
             "node_size" + std::to_string(node_size);
    }

    int k;                       ///< The number of nearest neighbors used to create the k-nearest
                                 ///< neighbor graph.
    int kc;                      ///< The scalar by which k is multiplied before querying the LSH
                                 ///< forest. The results are then ordered decreasing based on
                                 ///< linear-scan distances and the top k results returned.
    int fme_iterations;          ///< Maximum number of iterations of the fast multipole
                                 ///< embedder.
    bool fme_randomize;          ///< Whether or not to randomize the layout at the start.
    int fme_threads;             ///< The number of threads for the fast multipole embedder.
    int fme_precision;           ///< The number of coefficients of the multipole expansion.
    int sl_repeats;              ///< The number of repeats of the scaling layout algorithm.
    int sl_extra_scaling_steps;  ///< Sets the number of repeats of the scaling.
    double sl_scaling_min;       ///< The minimum scaling factor.
    double sl_scaling_max;       ///< The maximum scaling factor.
    ScalingType sl_scaling_type; ///< Defines the (relative) scale of the graph.
    int mmm_repeats;             ///< Number of repeats of the per-level layout algorithm.
    Placer placer;               ///< The  method  by  which  the  initial  positons  of  the
                                 ///< vertices  at  eachlevel are defined.
    Merger merger;               ///< The vertex merging strategy applied during the coarsening
                                 ///< phaseof the multilevel algorithm.
    double merger_factor;        ///< The ratio of the sizes between two levels up to
                                 ///< which the mergingis run.  Does not apply to all
                                 ///< merging strategies.
    int merger_adjustment;       ///< The  edge  length  adjustment  of  the  merging
                                 ///< algorithm.   Does  notapply to all merging
                                 ///< strategies.
    float node_size;             ///< The size of the nodes, which affects the magnitude of
                                 ///< their repellingforce. Decreasing  this  value  generally
                                 ///< resolves  overlaps  in  a  verycrowded tree.
  };

  /**
   * @brief The properties of a generated graph. An instance of this struct is
   * returned from the layout functions.
   *
   */
  struct GraphProperties
  {
    float mst_weight = 0.0;              ///< The total weight of the created spanning tree.
    uint32_t n_connected_components = 0; ///< The number of connected components.
    uint32_t n_isolated_vertices = 0;    ///< The number of isolated (lone) vertices.
    std::vector<uint32_t> degrees;       ///< The degrees of the vertices in the graph.
    std::vector<std::vector<std::pair<uint32_t, float>>>
        adjacency_list; ///< The adjacency list of the spanning tree.
    std::vector<std::vector<std::pair<uint32_t, float>>>
        adjacency_list_knn; ///< The adjacency list of the knn graph.
  };

  /**
   * @brief Genereates coordinates, edges and properties of a MST (via a kNN
   * graph) from an LSHForest instance.
   *
   * @param lsh_forest An LSHForest instance which is used to construct the kNN
   * graph.
   * @param config A LayoutConfiguration instance.
   * @param keep_knn Whether to keep the knn graph information as an adjacency list in GraphProperties.
   * @param create_mst Whether to create an MST before laying out the graph.
   * @param clear_lsh_forest Whether to clear the LSHForest after it's use (might
   * save memory).
   * @return std::tuple<std::vector<float>, std::vector<float>,
   * std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
   */
  std::tuple<std::vector<float>,
             std::vector<float>,
             std::vector<uint32_t>,
             std::vector<uint32_t>,
             GraphProperties>
  LayoutFromLSHForest(LSHForest &lsh_forest,
                      LayoutConfiguration config = LayoutConfiguration(),
                      bool keep_knn = false,
                      bool create_mst = true,
                      bool clear_lsh_forest = false);

  /**
   * @brief Generates an MST (via a kNN graph) from an LSHForest instance.
   *
   * @param lsh_forest An LSHForest instance which is used to construct the kNN
   * graph.
   * @param k The number of nearest neighbors used to create the kNN graph.
   * @param kc The factor by which k is multiplied when retrieving nearest
   * neighbors.
   * @return std::tuple<std::vector<uint32_t>, std::vector<uint32_t>, std::vector<float>>
   */
  std::tuple<std::vector<uint32_t>, std::vector<uint32_t>, std::vector<float>>
  MSTFromLSHForest(LSHForest &lsh_forest, uint32_t k, uint32_t kc = 10);

  /**
   * @brief Generates an MST from an edge list.
   *
   * @param vertex_count The number of vertices in the input graph.
   * @param edges An edge list in the form of [(from, to, weight)].
   * @return std::tuple<std::vector<uint32_t>, std::vector<uint32_t>, std::vector<float>>
   */
  std::tuple<std::vector<uint32_t>, std::vector<uint32_t>, std::vector<float>>
  MSTFromEdgeList(uint32_t vertex_count, const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges);

  /**
   * @brief Genereates coordinates, edges and properties of a MST from an edge
   * list.
   *
   * @param vertex_count The number of vertices in the input graph.
   * @param edges An edge list in the form of [(from, to, weight)].
   * @param config A LayoutConfiguration instance.
   * @param keep_knn Whether to keep the knn graph information as an adjacency list in GraphProperties.
   * @param create_mst Whether to create an MST before laying out the graph.
   * @return std::tuple<std::vector<float>, std::vector<float>,
   * std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
   */
  std::tuple<std::vector<float>,
             std::vector<float>,
             std::vector<uint32_t>,
             std::vector<uint32_t>,
             GraphProperties>
  LayoutFromEdgeList(
      uint32_t vertex_count,
      const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges,
      LayoutConfiguration config = LayoutConfiguration(),
      bool keep_knn = false,
      bool create_mst = true);

  /**
   * @brief Laying out an OGDF graph.
   *
   * @param g An OGDF Graph instance
   * @param vertex_count The number of vertices in the graph.
   * @param config A LayoutConfiguration instance.
   * @param gp An instance of a GraphProperties struct.
   * @return std::tuple<std::vector<float>, std::vector<float>,
   * std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
   */
  std::tuple<std::vector<float>,
             std::vector<float>,
             std::vector<uint32_t>,
             std::vector<uint32_t>,
             GraphProperties>
  LayoutInternal(ogdf::EdgeWeightedGraph<float> &g,
                 uint32_t vertex_count,
                 LayoutConfiguration config,
                 GraphProperties &gp);

  /**
   * @brief Calculates the quality of the vertex based on the actual nearest neighbors and the topological distances in the spanning tree.
   *
   * @param gp A GraphProperties object.
   * @param v The id of a vertex / data point.
   * @return std::vector<std::tuple<uint32_t, float, uint32_t>>
   */
  std::vector<std::tuple<uint32_t, float, uint32_t>>
  VertexQuality(GraphProperties &gp, uint32_t v);

  /**
   * @brief Calculates the mean quality of all vertices based on the actual nearest neighbors and the topological distances in the spanning tree.
   *
   * @param gp A GraphProperties object.
   * @return std::vector<float>
   */
  std::vector<float>
  MeanQuality(GraphProperties &gp);

  /**
   * @brief Gets the topological distances of a vertex to all other vertices.
   *
   * @param gp A GraphProperties object.
   * @param v The id of a vertex / data point.
   * @return std::vector<uint32_t>
   */
  std::vector<uint32_t>
  GetTopologicalDistances(GraphProperties &gp, uint32_t v);

  /**
   * @brief Creates an edge list from x, y coordinates and edge indices.
   *
   * @param x The x coordinates
   * @param y The y coordinates
   * @param s The indices of the from vertices
   * @param t The indices of the to vertices
   * @return std::tuple<std::vector<float>, std::vector<float>,
   * std::vector<float>, std::vector<float>>
   */
  std::tuple<std::vector<float>,
             std::vector<float>,
             std::vector<float>,
             std::vector<float>>
  MakeEdgeList(std::vector<float> x, std::vector<float> y,
               std::vector<uint32_t> s, std::vector<uint32_t> t);

}; // namespace tmap
#endif