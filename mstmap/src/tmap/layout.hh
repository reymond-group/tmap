/**
 * @file layout.hh
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Functions used for generating graph layouts from LSHForest instances and edge lists.
 * @version 0.1
 * @date 2019-06-17
 * 
 */

#ifndef LAYOUT_H
#define LAYOUT_H

#include <vector>
#include <tuple>
#include <stdint.h>

#include <ogdf/basic/Graph.h>
#include <ogdf/fileformats/GraphIO.h>
#include <ogdf/basic/GraphAttributes.h>
#include <ogdf/basic/extended_graph_alg.h>

#include <ogdf/energybased/FastMultipoleEmbedder.h>
#include <ogdf/energybased/multilevel_mixer/ScalingLayout.h>

#include <ogdf/energybased/multilevel_mixer/BarycenterPlacer.h>
#include <ogdf/energybased/multilevel_mixer/SolarPlacer.h>
#include <ogdf/energybased/multilevel_mixer/CirclePlacer.h>
#include <ogdf/energybased/multilevel_mixer/MedianPlacer.h>
#include <ogdf/energybased/multilevel_mixer/RandomPlacer.h>
#include <ogdf/energybased/multilevel_mixer/ZeroPlacer.h>

#include <ogdf/energybased/multilevel_mixer/SolarMerger.h>
#include <ogdf/energybased/multilevel_mixer/EdgeCoverMerger.h>
#include <ogdf/energybased/multilevel_mixer/IndependentSetMerger.h>
#include <ogdf/energybased/multilevel_mixer/LocalBiconnectedMerger.h>

#include <ogdf/basic/PreprocessorLayout.h>
#include <ogdf/packing/ComponentSplitterLayout.h>
#include <ogdf/packing/TileToRowsCCPacker.h>

#include <ogdf/energybased/FMMMLayout.h>
#include <ogdf/energybased/SpringEmbedderGridVariant.h>

#include <iostream>
#include <fstream>
#include <typeinfo>

// Forward declaration of LSHForest.
class LSHForest;

/**
 * @brief The placers available in OGDF.
 * 
 */
enum struct Placer
{
    Barycenter = 0,
    Solar = 1,
    Circle = 2,
    Median = 3,
    Random = 4,
    Zero = 5
};
static const std::string placer_values[] = {"Barycenter", "Solar", "Circle", "Median", "Random", "Zero"};

/**
 * @brief The mergers available in OGDF.
 * 
 */
enum struct Merger
{
    EdgeCover = 0,
    LocalBiconnected = 1,
    Solar = 2,
    IndependentSet = 3
};
static const std::string merger_values[] = {"EdgeCover", "LocalBiconnected", "Solar", "IndependentSet"};

/**
 * @brief The scaling types available in OGDF.
 * 
 */
enum struct ScalingType
{
    Absolute = 0,
    RelativeToAvgLength = 1,
    RelativeToDesiredLength = 2,
    RelativeToDrawing = 3
};
static const std::string scaling_types_values[] = {"Absolute", "RelativeToAvgLength", "RelativeToDesiredLength", "RelativeToDrawing"};

/**
 * @brief A struct containing all the configuration options available for and applied to a layout.
 * 
 */
struct LayoutConfiguration
{
    /**
     * @brief Construct a new Layout Configuration object.
     * 
     */
    LayoutConfiguration() : k(10), kc(10), fme_iterations(1000), fme_randomize(false), fme_threads(4), fme_precision(4),
                            sl_repeats(1), sl_extra_scaling_steps(1), sl_scaling_min(5.0), sl_scaling_max(20.0),
                            sl_scaling_type(ScalingType::RelativeToDrawing),
                            mmm_repeats(1),
                            placer(Placer::Barycenter),
                            merger(Merger::LocalBiconnected), merger_factor(2.0), merger_adjustment(0),
                            node_size(1.0) {}
    
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
               "sl_extra_scaling_steps: " + std::to_string(sl_extra_scaling_steps) + '\n' +
               "sl_scaling_x: " + std::to_string(sl_scaling_min) + '\n' +
               "sl_scaling_y: " + std::to_string(sl_scaling_max) + '\n' +
               "sl_scaling_type: " + scaling_types_values[(int)sl_scaling_type] + '\n' +
               "mmm_repeats: " + std::to_string(mmm_repeats) + '\n' +
               "placer: " + placer_values[(int)placer] + '\n' +
               "merger: " + merger_values[(int)merger] + '\n' +
               "merger_factor: " + std::to_string(merger_factor) + '\n' +
               "merger_adjustment: " + std::to_string(merger_adjustment)  + '\n' +
               "node_size" + std::to_string(node_size);
    }

    int k;
    int kc;
    int fme_iterations;
    bool fme_randomize;
    int fme_threads;
    int fme_precision;

    int sl_repeats;
    int sl_extra_scaling_steps;
    double sl_scaling_min;
    double sl_scaling_max;
    ScalingType sl_scaling_type;
    int mmm_repeats;
    Placer placer;
    Merger merger;
    double merger_factor;
    int merger_adjustment;
    float node_size;
};

/**
 * @brief The properties of a generated graph. An instance of this struct is returned from the layout functions.
 * 
 */
struct GraphProperties
{
    float mst_weight = 0.0;
    uint32_t n_connected_components = 0;
    uint32_t n_isolated_vertices = 0;
    std::vector<uint32_t> degrees;
    std::vector<std::vector<uint32_t>> adjacency_list;
};

/**
 * @brief Genereates coordinates, edges and properties of a MST (via a kNN graph) from an LSHForest instance.
 * 
 * @param lsh_forest An LSHForest instance which is used to construct the kNN graph.
 * @param config A LayoutConfiguration instance.
 * @param create_mst Whether to create an MST before laying out the graph.
 * @param clear_lsh_forest Whether to clear the LSHForest after it's use (might save memory).
 * @param weighted Whether the LSHForest instance contains weighted MinHash data.
 * @return std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties> 
 */
std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
LayoutFromLSHForest(LSHForest &lsh_forest, LayoutConfiguration config = LayoutConfiguration(), 
       bool create_mst = true, bool clear_lsh_forest = false, bool weighted = false);

/**
 * @brief Generates an MST (via a kNN graph) from an LSHForest instance.
 * 
 * @param lsh_forest An LSHForest instance which is used to construct the kNN graph.
 * @param k The number of nearest neighbors used to create the kNN graph.
 * @param kc The factor by which k is multiplied when retrieving nearest neighbors.
 * @param weighted Whether the LSHForest instance contains weighted MinHash data.
 * @return std::tuple<std::vector<uint32_t>, std::vector<uint32_t>> 
 */
std::tuple<std::vector<uint32_t>, std::vector<uint32_t>>
MSTFromLSHForest(LSHForest &lsh_forest, uint32_t k, uint32_t kc = 10, bool weighted = false);

/**
 * @brief Genereates coordinates, edges and properties of a MST from an edge list.
 * 
 * @param vertex_count The number of vertices in the input graph.
 * @param edges An edge list in the form of [(from, to, weight)].
 * @param config A LayoutConfiguration instance.
 * @param create_mst Whether to create an MST before laying out the graph.
 * @return std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties> 
 */
std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
LayoutFromEdgeList(uint32_t vertex_count, const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges,
       LayoutConfiguration config = LayoutConfiguration(), bool create_mst = true);

/**
 * @brief Laying out an OGDF graph.
 * 
 * @param g An OGDF Graph instance
 * @param vertex_count The number of vertices in the graph.
 * @param config A LayoutConfiguration instance.
 * @param gp An instance of a GraphProperties struct.
 * @return std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties> 
 */
std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>, GraphProperties>
LayoutInternal(ogdf::EdgeWeightedGraph<float> &g, uint32_t vertex_count, LayoutConfiguration config, GraphProperties &gp);

#endif