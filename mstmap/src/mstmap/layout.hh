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

#include <iostream>
#include <fstream>
#include <typeinfo>

// Forward declaration of LSHForest.
class LSHForest;

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
enum struct Merger
{
    EdgeCover = 0,
    LocalBiconnected = 1,
    Solar = 2,
    IndependentSet = 3
};
static const std::string merger_values[] = {"EdgeCover", "LocalBiconnected", "Solar", "IndependentSet"};
enum struct ScalingType
{
    Absolute = 0,
    RelativeToAvgLength = 1,
    RelativeToDesiredLength = 2,
    RelativeToDrawing = 3
};
static const std::string scaling_types_values[] = {"Absolute", "RelativeToAvgLength", "RelativeToDesiredLength", "RelativeToDrawing"};

struct LayoutConfiguration
{
    LayoutConfiguration() : k(10), kc(10), fme_iterations(1000), fme_randomize(false), fme_threads(4),
                            sl_repeats(1), sl_extra_scaling_steps(1), sl_scaling_x(5.0), sl_scaling_y(20.0),
                            sl_scaling_type(ScalingType::RelativeToDesiredLength),
                            mmm_repeats(1),
                            placer(Placer::Barycenter),
                            merger(Merger::LocalBiconnected), merger_factor(2.0), merger_adjustment(0),
                            node_size(1.0) {}
    std::string ToString() const
    {
        return std::string("k: ") + std::to_string(k) + '\n' +
               "kc: " + std::to_string(kc) + '\n' +
               "fme_iterations: " + std::to_string(fme_iterations) + '\n' +
               "fme_randomize: " + std::to_string(fme_randomize) + '\n' +
               "fme_threads: " + std::to_string(fme_threads) + '\n' +
               "sl_repeats: " + std::to_string(sl_repeats) + '\n' +
               "sl_extra_scaling_steps: " + std::to_string(sl_extra_scaling_steps) + '\n' +
               "sl_scaling_x: " + std::to_string(sl_scaling_x) + '\n' +
               "sl_scaling_y: " + std::to_string(sl_scaling_y) + '\n' +
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
    int sl_repeats;
    int sl_extra_scaling_steps;
    double sl_scaling_x;
    double sl_scaling_y;
    ScalingType sl_scaling_type;
    int mmm_repeats;
    Placer placer;
    Merger merger;
    double merger_factor;
    int merger_adjustment;
    float node_size;
};

std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>>
LayoutFromLSHForest(LSHForest &lsh_forest, LayoutConfiguration config = LayoutConfiguration(), 
       bool create_mst = true, bool clear_lsh_forest = false, bool weighted = false);

std::tuple<std::vector<uint32_t>, std::vector<uint32_t>>
MSTFromLSHForest(LSHForest &lsh_forest, uint32_t k, uint32_t kc = 10, bool weighted = false);

std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>>
LayoutFromEdgeList(uint32_t vertex_count, const std::vector<std::tuple<uint32_t, uint32_t, float>> &edges,
       LayoutConfiguration config = LayoutConfiguration(), bool create_mst = true);

std::tuple<std::vector<float>, std::vector<float>, std::vector<uint32_t>, std::vector<uint32_t>>
LayoutInternal(ogdf::Graph &g, uint32_t vertex_count, LayoutConfiguration config);

#endif