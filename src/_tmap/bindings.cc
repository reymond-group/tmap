/**
 * @file bindings.cc
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Pybind11 bindings for tmap.
 * @version 0.1
 * @date 2019-06-17
 *
 */

#include "layout.hh"
#include "lshforest.hh"
#include "minhash.hh"
#include "analyse.hh"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <stdint.h>
#include <tuple>
#include <vector>

using namespace tmap;

namespace py = pybind11;

// Used to convert python lists to STL vectors
template <class T>
std::vector<T> convert_list(py::list &l)
{
    std::vector<T> vec(l.size());
    for (size_t i = 0; i < l.size(); i++)
        vec[i] = l[i].cast<T>();
    return vec;
}

// Used to convert python lists of lists to STL vectors of vectors
template <class T>
std::vector<std::vector<T>> convert_list_of_lists(py::list &l)
{
    std::vector<std::vector<T>> vecs(l.size());
    for (size_t i = 0; i < l.size(); i++)
    {
        py::list sl = (py::list)l[i];
        vecs[i] = std::vector<T>(sl.size());
        for (size_t j = 0; j < sl.size(); j++)
            vecs[i][j] = sl[j].cast<T>();
    }

    return vecs;
}

// Used to convert python 2D arrays to STL vectors
template <class T>
std::vector<std::vector<T>> convert_array(py::array_t<T> arr)
{
    py::buffer_info buffi = arr.request();
    T *ptr = (T *)buffi.ptr;
    int rows = buffi.shape[0];
    int cols = buffi.shape[1];

    std::vector<std::vector<T>> vecs(rows);

    for (size_t i = 0; i < rows; i++)
    {
        vecs[i] = std::vector<T>(cols);
        for (size_t j = 0; j < cols; j++)
        {
            vecs[i][j] = ptr[i * cols + j];
        }
    }

    return vecs;
}

// Used to convert python 2D arrays to STL vectors, including a type change
template <class T>
std::vector<std::vector<float>> convert_array_to_float(py::array_t<T> arr)
{
    py::buffer_info buffi = arr.request();
    T *ptr = (T *)buffi.ptr;
    int rows = buffi.shape[0];
    int cols = buffi.shape[1];

    // std::cout << rows << std::endl;
    // std::cout << cols << std::endl;

    std::vector<std::vector<float>> vecs(rows);

    for (size_t i = 0; i < rows; i++)
    {
        vecs[i] = std::vector<float>(cols);
        for (size_t j = 0; j < cols; j++)
        {
            vecs[i][j] = (float)ptr[j * rows + i];
        }

        // std::cout << i << std::endl;
    }

    // std::cout << "done" << std::endl;
    // std::cout << rows << std::endl;
    // std::cout << cols << std::endl;

    return vecs;
}

class TestSuper
{
public:
    int x;
};

class TestSub : public TestSuper
{
public:
    int y;
};

// Extending the Minhash class to allow for passing lists as well as
// opaque types
class PyMinhash : public Minhash
{
public:
    using Minhash::BatchFromBinaryArray;
    using Minhash::BatchFromSparseBinaryArray;
    using Minhash::BatchFromWeightArray;
    using Minhash::FromBinaryArray;
    using Minhash::FromSparseBinaryArray;
    using Minhash::FromWeightArray;
    using Minhash::Minhash;

    std::vector<uint32_t> FromBinaryArray(py::list &list)
    {
        std::vector<uint8_t> vec = convert_list<uint8_t>(list);
        return Minhash::FromBinaryArray(vec);
    }

    std::vector<std::vector<uint32_t>> BatchFromBinaryArray(py::list &list)
    {
        std::vector<std::vector<uint8_t>> vecs = convert_list_of_lists<uint8_t>(list);
        return Minhash::BatchFromBinaryArray(vecs);
    }

    std::vector<std::vector<uint32_t>> BatchFromBinaryArray(py::array_t<uint8_t> &arr)
    {
        std::vector<std::vector<uint8_t>> vecs = convert_array<uint8_t>(arr);
        return Minhash::BatchFromBinaryArray(vecs);
    }

    std::vector<uint32_t> FromSparseBinaryArray(py::list &list)
    {
        std::vector<uint32_t> vec = convert_list<uint32_t>(list);
        return Minhash::FromSparseBinaryArray(vec);
    }

    std::vector<std::vector<uint32_t>> BatchFromSparseBinaryArray(py::list &list)
    {
        std::vector<std::vector<uint32_t>> vecs = convert_list_of_lists<uint32_t>(list);
        return Minhash::BatchFromSparseBinaryArray(vecs);
    }

    std::vector<std::vector<uint32_t>> BatchFromSparseBinaryArray(py::array_t<uint32_t> &arr)
    {
        std::vector<std::vector<uint32_t>> vecs = convert_array<uint32_t>(arr);
        return Minhash::BatchFromSparseBinaryArray(vecs);
    }

    std::vector<uint32_t> FromWeightArray(py::list &list, const std::string &method = "ICWS")
    {
        std::vector<float> vec = convert_list<float>(list);
        return Minhash::FromWeightArray(vec, method);
    }

    std::vector<std::vector<uint32_t>> BatchFromWeightArray(py::list &list, const std::string &method = "ICWS")
    {
        std::vector<std::vector<float>> vecs = convert_list_of_lists<float>(list);
        return Minhash::BatchFromWeightArray(vecs, method);
    }

    std::vector<std::vector<uint32_t>> BatchFromWeightArray(py::array_t<float, 32> &arr, const std::string &method = "ICWS")
    {
        std::vector<std::vector<float>> vecs = convert_array<float>(arr);
        return Minhash::BatchFromWeightArray(vecs, method);
    }

    std::vector<std::vector<uint32_t>> BatchFromWeightArray(py::array_t<double> &arr, const std::string &method = "ICWS")
    {
        std::vector<std::vector<float>> vecs = convert_array_to_float<double>(arr);
        return Minhash::BatchFromWeightArray(vecs, method);
    }
};

// Defining opaque types (see bindings below as well)
PYBIND11_MAKE_OPAQUE(std::vector<uint8_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint16_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint32_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint64_t>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);

// In order to make the python module useable from R via reticulate,
// layout methods that return native python lists rather than bound
// std vectors are implemented.
py::tuple LayoutFromLSHForestNative(LSHForest &lsh_forest,
                                    LayoutConfiguration config = LayoutConfiguration(),
                                    bool keep_knn = false,
                                    bool create_mst = true,
                                    bool clear_lsh_forest = false)
{
    auto result = LayoutFromLSHForest(lsh_forest, config, keep_knn, create_mst, clear_lsh_forest);
    py::list x = py::cast(std::get<0>(result));
    py::list y = py::cast(std::get<1>(result));
    py::list s = py::cast(std::get<2>(result));
    py::list t = py::cast(std::get<3>(result));
    py::object gp = py::cast(std::get<4>(result));
    return py::make_tuple(x, y, s, t, gp);
}

py::tuple LayoutFromEdgeListNative(uint32_t vertex_count,
                                   std::vector<std::tuple<uint32_t, uint32_t, float>> &edges,
                                   LayoutConfiguration config = LayoutConfiguration(),
                                   bool keep_knn = false,
                                   bool create_mst = true)
{
    auto result = LayoutFromEdgeList(vertex_count, edges, config, keep_knn, create_mst);
    py::list x = py::cast(std::get<0>(result));
    py::list y = py::cast(std::get<1>(result));
    py::list s = py::cast(std::get<2>(result));
    py::list t = py::cast(std::get<3>(result));
    py::object gp = py::cast(std::get<4>(result));
    return py::make_tuple(x, y, s, t, gp);
}

py::tuple MakeEdgeListNative(std::vector<float> x, std::vector<float> y,
                             std::vector<uint32_t> s, std::vector<uint32_t> t)
{
    auto result = MakeEdgeList(x, y, s, t);
    py::list x1 = py::cast(std::get<0>(result));
    py::list y1 = py::cast(std::get<1>(result));
    py::list x2 = py::cast(std::get<2>(result));
    py::list y2 = py::cast(std::get<3>(result));
    return py::make_tuple(x1, y1, x2, y2, py::cast(x), py::cast(y));
}

// Windows and Mac compiles complained about dtype being non const

// Make life easier for R people
template <class T>
py::tuple map(py::array_t<T> arr, uint32_t dims = 128, uint32_t n_trees = 8,
              const std::string &dtype = "binary",
              LayoutConfiguration config = LayoutConfiguration(),
              bool file_backed = false,
              unsigned int seed = 42)
{
    py::tuple result;
    if (dtype == "binary")
    {
        Minhash mh(dims);
        LSHForest lf(dims, n_trees, true, file_backed, false);

        auto vecs = convert_array<uint8_t>(arr);
        auto hashes = mh.BatchFromBinaryArray(vecs);
        lf.BatchAdd(hashes);
        lf.Index();
        result = LayoutFromLSHForestNative(lf, config);
    }
    else if (dtype == "sparse")
    {
        Minhash mh(dims);
        LSHForest lf(dims, n_trees, true, file_backed, false);

        auto vecs = convert_array<uint32_t>(arr);
        auto hashes = mh.BatchFromSparseBinaryArray(vecs);
        lf.BatchAdd(hashes);
        lf.Index();
        result = LayoutFromLSHForestNative(lf, config);
    }
    else if (dtype == "weighted")
    {
        auto vecs = convert_array_to_float<T>(arr);

        Minhash mh(vecs[0].size(), seed, dims);
        LSHForest lf(dims * 2, n_trees, true, file_backed, true);

        auto hashes = mh.BatchFromWeightArray(vecs);
        std::cout << "Have hashes" << std::endl;
        std::cout << vecs[0].size() << std::endl;
        lf.BatchAdd(hashes);
        std::cout << "added to lf" << std::endl;
        lf.Index();
        std::cout << "indexed" << std::endl;
        auto r = LayoutFromLSHForest(lf);
        std::cout << "have result in" << std::endl;
        result = LayoutFromLSHForestNative(lf, config);
        std::cout << "have result" << std::endl;
    }
    else
    {
        throw std::invalid_argument("dtype has to be 'binary', 'sparse', or 'weighted'");
    }

    return result;
}

PYBIND11_MODULE(_tmap, m)
{
    py::bind_vector<std::vector<uint8_t>>(
        m, "VectorUchar", "Unsigned 8-bit int vector.");
    py::bind_vector<std::vector<uint16_t>>(
        m, "VectorUsmall", "Unsigned 16-bit int vector.");
    py::bind_vector<std::vector<uint32_t>>(
        m, "VectorUint", "Unsigned 32-bit int vector.");
    py::bind_vector<std::vector<float>>(
        m, "VectorFloat", "Unsigned 32-bit float vector.");
    py::bind_vector<std::vector<uint64_t>>(
        m, "VectorUlong", "Unsigned 64-bit int vector.");

    py::enum_<ScalingType>(m, "ScalingType", py::arithmetic(), R"pbdoc(
        The scaling types available in OGDF. The class is to be used as an enum.

        Notes:
            The available values are

            :obj:`ScalingType.Absolute`: Absolute factor, can be used to scale relative to level size change.
            
            :obj:`ScalingType.RelativeToAvgLength`: Scales by a factor relative to the average edge weights.
            
            :obj:`ScalingType.RelativeToDesiredLength`: Scales by a factor relative to the disired edge length.
            
            :obj:`ScalingType.RelativeToDrawing`: Scales by a factor relative to the drawing.
    )pbdoc")
        .value("Absolute", ScalingType::Absolute)
        .value("RelativeToAvgLength", ScalingType::RelativeToAvgLength)
        .value("RelativeToDesiredLength", ScalingType::RelativeToDesiredLength)
        .value("RelativeToDrawing", ScalingType::RelativeToDrawing)
        .export_values();

    py::enum_<Placer>(m, "Placer", py::arithmetic(), R"pbdoc(
        The places available in OGDF. The class is to be used as an enum.

        Notes:
            The available values are

            :obj:`Placer.Barycenter`: Places a vertex at the barycenter of its neighbors' position.
            
            :obj:`Placer.Solar`: Uses information of the merging phase of the solar merger. Places a new vertex on the direct line between two suns.
            
            :obj:`Placer.Circle`: Places the vertices in a circle around the barycenter and outside of the current drawing
            
            :obj:`Placer.Median`: Places a vertex at the median position of the neighbor nodes for each coordinate axis.
            
            :obj:`Placer.Random`: Places a vertex at a random position within the smallest circle containing all vertices around the barycenter of the current drawing.
            
            :obj:`Placer.Zero`: Places a vertex at the same position as its representative in the previous level.
    )pbdoc")
        .value("Barycenter", Placer::Barycenter)
        .value("Solar", Placer::Solar)
        .value("Circle", Placer::Circle)
        .value("Median", Placer::Median)
        .value("Random", Placer::Random)
        .value("Zero", Placer::Zero)
        .export_values();

    py::enum_<Merger>(m, "Merger", py::arithmetic(), R"pbdoc(
        The mergers available in OGDF. The class is to be used as an enum.

        Notes:
            The available values are

            :obj:`Merger.EdgeCover`: Based on the matching merger. Computes an edge cover such that each contained edge is incident to at least one unmatched vertex. The cover edges are then used to merge their adjacent vertices.
            
            :obj:`Merger.LocalBiconnected`: Based on the edge cover merger. Avoids distortions by checking whether biconnectivity will be lost in the local neighborhood around the potential merging position.
            
            :obj:`Merger.Solar`: Vertices are partitioned into solar systems, consisting of sun, planets and moons. The systems are then merged into the sun vertices.
            
            :obj:`Merger.IndependentSet`: Uses a maximal independent set filtration. See GRIP for details.
    )pbdoc")
        .value("EdgeCover", Merger::EdgeCover)
        .value("LocalBiconnected", Merger::LocalBiconnected)
        .value("Solar", Merger::Solar)
        .value("IndependentSet", Merger::IndependentSet)
        .export_values();

    py::class_<LayoutConfiguration>(m, "LayoutConfiguration", R"pbdoc(
        A container for configuration options for :obj:`layout_from_lsh_forest()` and :obj:`layout_from_edge_list()`.

        Attributes:
            int k (:obj:`int`): The number of nearest neighbors used to create the k-nearest neighbor graph.
            int kc (:obj:`int`): The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned.
            int fme_iterations (:obj:`int`): Maximum number of iterations of the fast multipole embedder.
            bool fme_randomize (:obj:`bool`): Whether or not to randomize the layout at the start.
            int fme_threads (:obj:`int`): The number of threads for the fast multipole embedder.
            int fme_precision (:obj:`int`): The number of coefficients of the multipole expansion.
            int sl_repeats (:obj:`int`): The number of repeats of the scaling layout algorithm.
            int sl_extra_scaling_steps (:obj:`int`): Sets the number of repeats of the scaling.
            double sl_scaling_min (:obj:`float`): The minimum scaling factor.
            double sl_scaling_max (:obj:`float`): The maximum scaling factor.
            ScalingType sl_scaling_type (:obj:`ScalingType`): Defines the (relative) scale of the graph.
            int mmm_repeats (:obj:`int`): Number of repeats of the per-level layout algorithm.
            Placer placer (:obj:`Placer`): The  method  by  which  the  initial  positons  of  the  vertices  at  eachlevel are defined.
            Merger merger (:obj:`Merger`): The vertex merging strategy applied during the coarsening phaseof the multilevel algorithm.
            double merger_factor (:obj:`float`): The ratio of the sizes between two levels up to which the mergingis run.  Does not apply to all merging strategies.
            int merger_adjustment (:obj:`int`): The  edge  length  adjustment  of  the  merging  algorithm.   Does  notapply to all merging strategies.
            float node_size (:obj:`float`): The size of the nodes, which affects the magnitude of their repellingforce. Decreasing  this  value  generally  resolves  overlaps  in  a  verycrowded tree.
    )pbdoc")
        .def(py::init(), R"pbdoc(
            Constructor for the class :obj:`LayoutConfiguration`.
        )pbdoc")
        .def_readwrite("k", &LayoutConfiguration::k)
        .def_readwrite("kc", &LayoutConfiguration::kc)
        .def_readwrite("fme_iterations", &LayoutConfiguration::fme_iterations)
        .def_readwrite("fme_randomize", &LayoutConfiguration::fme_randomize)
        .def_readwrite("fme_threads", &LayoutConfiguration::fme_threads)
        .def_readwrite("fme_precision", &LayoutConfiguration::fme_precision)
        .def_readwrite("sl_repeats", &LayoutConfiguration::sl_repeats)
        .def_readwrite("sl_extra_scaling_steps",
                       &LayoutConfiguration::sl_extra_scaling_steps)
        .def_readwrite("sl_scaling_min", &LayoutConfiguration::sl_scaling_min)
        .def_readwrite("sl_scaling_max", &LayoutConfiguration::sl_scaling_max)
        .def_readwrite("sl_scaling_type", &LayoutConfiguration::sl_scaling_type)
        .def_readwrite("mmm_repeats", &LayoutConfiguration::mmm_repeats)
        .def_readwrite("placer", &LayoutConfiguration::placer)
        .def_readwrite("merger", &LayoutConfiguration::merger)
        .def_readwrite("merger_factor", &LayoutConfiguration::merger_factor)
        .def_readwrite("merger_adjustment", &LayoutConfiguration::merger_adjustment)
        .def_readwrite("node_size", &LayoutConfiguration::node_size)
        .def("__repr__", &LayoutConfiguration::ToString);

    py::class_<GraphProperties>(m, "GraphProperties", R"pbdoc(
        Contains properties of the minimum spanning tree (or forest) generated by :obj:`layout_from_lsh_forest()` and :obj:`layout_from_edge_list()`.

        Attributes:
            mst_weight (:obj:`float`): The total weight of the minimum spanning tree.
            n_connected_components (:obj:`int`): The number of connected components in the minimum spanning forest.
            n_isolated_vertices (:obj:`int`) The number of isolated vertices in the minimum spanning forest.
            degrees (:obj:`VectorUint`): The degrees of all vertices in the minimum spanning tree (or forest).
            adjacency_list(:obj:`List` of :obj:`VectorUint`): The adjaceny lists for all vertices in the minimum spanning tree (or forest).
    )pbdoc")
        .def(py::init(), R"pbdoc(
            Constructor for the class :obj:`GraphProperties`.
        )pbdoc")
        .def_readonly("mst_weight", &GraphProperties::mst_weight)
        .def_readonly("n_connected_components",
                      &GraphProperties::n_connected_components)
        .def_readonly("n_isolated_vertices", &GraphProperties::n_isolated_vertices)
        .def_readonly("degrees", &GraphProperties::degrees)
        .def_readonly("adjacency_list", &GraphProperties::adjacency_list)
        .def_readonly("adjacency_list_knn", &GraphProperties::adjacency_list_knn);

    m.def("map",
          &map<double>,
          py::arg("arr"),
          py::arg("dims") = 128,
          py::arg("n_trees") = 8,
          py::arg("dtype") = "binary",
          py::arg("config") = LayoutConfiguration(),
          py::arg("file_backed") = false,
          py::arg("seed") = 42,
          R"pbdoc(
        Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an :obj:`LSHForest` instance. This method returns native python lists and objects.
        
        Arguments:
            arr (:obj:`Array`): A numpy :obj:`Array` instance
        
        Keyword Arguments:
            dims (:obj:`int`, optional): The number of permutations to use for the MinHash algorithm
            n_trees (:obj:`int`, optional): The number of forests to use in the LSHForest data structure
            dtype (:obj:`str`, optional): The type of data that is supplied, can be 'binary', 'sparse', or 'weighted'
            config (:obj:`LayoutConfiguration`, optional): An :obj:`LayoutConfiguration` instance
            file_backed (:obj:`bool`) Whether to store the data on disk rather than in main memory (experimental)
            seed (:obj:`int`): The seed used for the random number generator(s)

        Returns:
            :obj:`Tuple[List, List, List, List, Object]` The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph
    )pbdoc");

    m.def("layout_from_lsh_forest",
          &LayoutFromLSHForest,
          py::arg("lsh_forest"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("keep_knn") = false,
          py::arg("create_mst") = true,
          py::arg("clear_lsh_forest") = false,
          R"pbdoc(
        Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an :obj:`LSHForest` instance.
        
        Arguments:
            lsh_forest (:obj:`LSHForest`): An :obj:`LSHForest` instance
        
        Keyword Arguments:
            config (:obj:`LayoutConfiguration`, optional): An :obj:`LayoutConfiguration` instance
            create_mst (:obj:`bool`, optional): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph
            clear_lsh_forest (:obj:`bool`, optional): Whether to run :obj:`clear()` on the :obj:`LSHForest` instance after k-nearest negihbor graph and MST creation and before layout

        Returns:
            :obj:`Tuple[VectorFloat, VectorFloat, VectorUint, VectorUint, GraphProperties]` The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph
    )pbdoc");

    m.def("layout_from_lsh_forest_native",
          &LayoutFromLSHForestNative,
          py::arg("lsh_forest"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("keep_knn") = false,
          py::arg("create_mst") = true,
          py::arg("clear_lsh_forest") = false,
          R"pbdoc(
        Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an :obj:`LSHForest` instance. This method returns native python lists and objects.
        
        Arguments:
            lsh_forest (:obj:`LSHForest`): An :obj:`LSHForest` instance
        
        Keyword Arguments:
            config (:obj:`LayoutConfiguration`, optional): An :obj:`LayoutConfiguration` instance
            create_mst (:obj:`bool`, optional): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph
            clear_lsh_forest (:obj:`bool`, optional): Whether to run :obj:`clear()` on the :obj:`LSHForest` instance after k-nearest negihbor graph and MST creation and before layout

        Returns:
            :obj:`Tuple[List, List, List, List, Object]` The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph
    )pbdoc");

    m.def("mst_from_lsh_forest",
          &MSTFromLSHForest,
          py::arg("lsh_forest"),
          py::arg("k"),
          py::arg("kc") = 10,
          R"pbdoc(
        Create minimum spanning tree topology from an :obj:`LSHForest` instance.
        
        Arguments:
            lsh_forest (:obj:`LSHForest`): An :obj:`LSHForest` instance
            int k (:obj:`int`): The number of nearest neighbors used to create the k-nearest neighbor graph
        
        Keyword Arguments:
            int kc (:obj:`int`, optional): The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned

        Returns:
            :obj:`Tuple[VectorUint, VectorUint, VectorFloat]`: the topology of the minimum spanning tree of the data indexed in the LSH forest
    )pbdoc");

    m.def("mst_from_edge_list",
          &MSTFromEdgeList,
          py::arg("vertex_count"),
          py::arg("edges"),
          R"pbdoc(
        Create minimum spanning tree topology from an edge list.
        
        Arguments:
            vertex_count (:obj:`int`): The number of vertices in the edge list
            edges (:obj:`List` of :obj:`Tuple[int, int, float]`): An edge list defining a graph

        Returns:
            :obj:`Tuple[VectorUint, VectorUint, VectorFloat]`: the topology of the minimum spanning tree of the data from the edge list
    )pbdoc");

    m.def("layout_from_edge_list",
          &LayoutFromEdgeList,
          py::arg("vertex_count"),
          py::arg("edges"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("keep_knn") = false,
          py::arg("create_mst") = true,
          R"pbdoc(
        Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an edge list.
        
        Arguments:
            vertex_count (:obj:`int`): The number of vertices in the edge list
            edges (:obj:`List` of :obj:`Tuple[int, int, float]`): An edge list defining a graph
        
        Keyword Arguments:
            config (:obj:`LayoutConfiguration`, optional): An :obj:`LayoutConfiguration` instance
            create_mst (:obj:`bool`): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph

        Returns:
            :obj:`Tuple[VectorFloat, VectorFloat, VectorUint, VectorUint, GraphProperties]`: The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph
    )pbdoc");

    m.def("layout_from_edge_list_native",
          &LayoutFromEdgeListNative,
          py::arg("vertex_count"),
          py::arg("edges"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("keep_knn") = false,
          py::arg("create_mst") = true,
          R"pbdoc(
        Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an edge list. This method returns native python lists and objects.
        
        Arguments:
            vertex_count (:obj:`int`): The number of vertices in the edge list
            edges (:obj:`List` of :obj:`Tuple[int, int, float]`): An edge list defining a graph
        
        Keyword Arguments:
            config (:obj:`LayoutConfiguration`, optional): An :obj:`LayoutConfiguration` instance
            create_mst (:obj:`bool`): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph

        Returns:
            :obj:`Tuple[VectorFloat, VectorFloat, VectorUint, VectorUint, GraphProperties]`: The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph
    )pbdoc");

    m.def("make_edge_list",
          &MakeEdgeList,
          py::arg("x"),
          py::arg("y"),
          py::arg("s"),
          py::arg("t"),
          R"pbdoc(
        Creates an edge list from x, y coordinates and edge indices.
        
        Arguments:
            x (:obj:`VectorFloat`): The x coordinates
            y (:obj:`VectorFloat`): The y coordinates
            s (:obj:`VectorUint`): The indices of the from vertices
            t (:obj:`VectorUint`): The indices of the to vertices
        
        Returns:
            :obj:`Tuple[VectorFloat, VectorFloat, VectorFloat, VectorFloat]`: Coordinates in edge list form
    )pbdoc");

    m.def("make_edge_list_native",
          &MakeEdgeListNative,
          py::arg("x"),
          py::arg("y"),
          py::arg("s"),
          py::arg("t"),
          R"pbdoc(
        brief Creates an edge list from x, y coordinates and edge indices. This method returns native python lists and objects. Also returns coordinates of vertices
        
        Arguments:
            x (:obj:`VectorFloat`): The x coordinates
            y (:obj:`VectorFloat`): The y coordinates
            s (:obj:`VectorUint`): The indices of the from vertices
            t (:obj:`VectorUint`): The indices of the to vertices
        
        Returns:
            :obj:`Tuple[List, List, List, List, List, List]`: Coordinates in edge list form and the vertex coordinates
    )pbdoc");

    m.def("vertex_quality",
          &VertexQuality,
          py::arg("gp"),
          py::arg("v"),
          R"pbdoc(
        brief Computes the visualization quality of a vertex based on it's true nearest neighbors and their distribution in the tree.
        
        Arguments:
            gp (:obj:`VectorFloat`): A GraphProperties object
            v (:obj:`int`): The vertex id of the node to be analysed
        
        Returns:
            :obj:`List[Tuple[int, float, int]]`: The qualities based on the degrees in the knn graph.
    )pbdoc");

    m.def("mean_quality",
          &MeanQuality,
          py::arg("gp"),
          R"pbdoc(
        brief Calculates the mean quality of all vertices based on the actual nearest neighbors and the topological distances in the spanning tree.
        
        Arguments:
            gp (:obj:`VectorFloat`): A GraphProperties object
        
        Returns:
            :obj:`List[float]`: The average topological distances ordered by k-nearest neighbor.
    )pbdoc");

    m.def("get_topological_distances",
          &GetTopologicalDistances,
          py::arg("gp"),
          py::arg("v"),
          R"pbdoc(
        brief Gets the topological distances of a vertex to all other vertices.
        
        Arguments:
            gp (:obj:`VectorFloat`): A GraphProperties object
            v (:obj:`int`): The vertex id of the node to be analysed
        
        Returns:
            :obj:`List[int]`: The topological distances to vertex v.
    )pbdoc");

    py::class_<LSHForest>(m, "LSHForest", R"pbdoc(
        A LSH forest data structure which incorporates optional linear scan to increase the recovery performance. Most query methods are available in parallelized versions named with a :obj:`batch_` prefix.
    )pbdoc")
        .def(py::init<unsigned int, unsigned int, bool, bool, bool>(),
             py::arg("d") = 128,
             py::arg("l") = 8,
             py::arg("store") = true,
             py::arg("file_backed") = false,
             py::arg("weighted") = false,
             R"pbdoc(
            Constructor for the class :obj:`LSHForest`.

            Keyword Arguments:
                d (:obj:`int`): The dimensionality of the MinHashe vectors to be added
                l (:obj:`int`): The number of prefix trees used when indexing data
                store (:obj:`bool`) Whether to store the added MinHash vectors. This is required when using linear scan in queries
                file_backed (:obj:`bool`) Whether to store the data on disk rather than in main memory (experimental)
        )pbdoc")
        .def("add", &LSHForest::Add, R"pbdoc(
            Add a MinHash vector to the LSH forest.

            Arguments:
                vecs (:obj:`VectorUint`): A MinHash vector that is to be added to the LSH forest
        )pbdoc")
        .def("batch_add", &LSHForest::BatchAdd, R"pbdoc(
            Add a list MinHash vectors to the LSH forest (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`VectorUint`): A list of MinHash vectors that is to be added to the LSH forest
        )pbdoc")
        .def("fit", &LSHForest::Fit, R"pbdoc(
            Add Minhashes with labels to this LSHForest (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`VectorUint`): A list of MinHash vectors that is to be added to the LSH forest
                labels (:obj:`VectorUint`) A vector containing labels.
        )pbdoc")
        .def("predict",
             &LSHForest::Predict,
             py::arg("vecs"),
             py::arg("k") = 10,
             py::arg("kc") = 10,
             py::arg("weighted") = false,
             R"pbdoc(
            Predict labels of Minhashes using the kNN algorithm (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`VectorUint`): A list of MinHash vectors that is to be added to the LSH forest
                k (:obj:`int`) The degree of the kNN algorithm
                kc (:obj:`int`) The scalar by which k is multiplied before querying the LSH
                weighted (:obj:`bool` Whether distances are used as weights by the knn algorithm)
            Returns:
                :obj:`VectorUint` The predicted labels
        )pbdoc")
        .def("index", &LSHForest::Index, R"pbdoc(
            Index the LSH forest. This has to be run after each time new MinHashes were added.
        )pbdoc")
        .def("is_clean", &LSHForest::IsClean, R"pbdoc(
            Returns a boolean indicating whether or not the LSH forest has been indexed after the last MinHash vector was added.

            Returns:
                :obj:`bool`: :obj:`True` if :obj:`index()` has been run since MinHash vectors have last been added using :obj:`add()` or :obj:`batch_add()`. :obj:`False` otherwise
        )pbdoc")
        .def("query_linear_scan",
             &LSHForest::QueryLinearScan,
             py::arg("vec"),
             py::arg("k"),
             py::arg("kc") = 10,
             R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query_linear_scan_exclude",
             &LSHForest::QueryLinearScanExclude,
             py::arg("vec"),
             py::arg("k"),
             py::arg("exclude"),
             py::arg("kc") = 10,
             R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                exclude (:obj:`List` of :obj:`VectorUint`) A list of ids of indexed MinHash vectors to be excluded from the search
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query_linear_scan_by_id",
             &LSHForest::QueryLinearScanById,
             py::arg("id"),
             py::arg("k"),
             py::arg("kc") = 10,
             R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                id (:obj:`int`): The id of an indexed MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query_linear_scan_exclude_by_id",
             &LSHForest::QueryLinearScanExcludeById,
             py::arg("id"),
             py::arg("k"),
             py::arg("exclude"),
             py::arg("kc") = 10,
             R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                id (:obj:`int`): The id of an indexed MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                exclude (:obj:`List` of :obj:`VectorUint`) A list of ids of indexed MinHash vectors to be excluded from the search
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("linear_scan",
             &LSHForest::LinearScan,
             py::arg("vec"),
             py::arg("indices"),
             py::arg("k") = 10,
             R"pbdoc(
            Query a subset of indexed MinHash vectors using linear scan.

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                indices (:obj:`VectorUint`) The ids of indexed MinHash vectors that define the subset to be queried

            Keyword arguments:
                k (:obj:`int`): The number of nearest neighbors to be retrieved

            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query", &LSHForest::Query, R"pbdoc(
            Query the LSH forest for k-nearest neighbors.

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Returns:
                :obj:`VectorUint`: The results of the query
        )pbdoc")
        .def("query_exclude", &LSHForest::QueryExclude, R"pbdoc(
            Query the LSH forest for k-nearest neighbors.

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                exclude (:obj:`List` of :obj:`VectorUint`) A list of ids of indexed MinHash vectors to be excluded from the search
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Returns:
                :obj:`VectorUint`: The results of the query
        )pbdoc")
        .def("query_by_id", &LSHForest::QueryById, R"pbdoc(
            Query the LSH forest for k-nearest neighbors.

            Arguments:
                id (:obj:`int`): The id of an indexed MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Returns:
                :obj:`VectorUint`: The results of the query
        )pbdoc")
        .def("query_exclude_by_id", &LSHForest::QueryExcludeById, R"pbdoc(
            Query the LSH forest for k-nearest neighbors.

            Arguments:
                id (:obj:`int`): The id of an indexed MinHash vector
                exclude (:obj:`List` of :obj:`VectorUint`) A list of ids of indexed MinHash vectors to be excluded from the search
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Returns:
                :obj:`VectorUint`: The results of the query
        )pbdoc")
        .def("batch_query", &LSHForest::BatchQuery, R"pbdoc(
            Query the LSH forest for k-nearest neighbors (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`VectorUint`): The query MinHash vectors
                k (:obj:`int`): The number of nearest neighbors to be retrieved

            Returns:
                :obj:`List` of :obj:`VectorUint`: The results of the queries
        )pbdoc")
        .def("get_all_nearest_neighbors",
             &LSHForest::GetAllNearestNeighbors,
             py::arg("k"),
             py::arg("kc") = 10,
             R"pbdoc(
            Get the k-nearest neighbors of all indexed MinHash vectors.

            Arguments:
                k (:obj:`int`): The number of nearest neighbors to be retrieved

            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival

            Returns:
                :obj:`VectorUint` The ids of all k-nearest neighbors
        )pbdoc")
        .def("get_knn_graph",
             &LSHForest::GetKNNGraph,
             py::arg("from"),
             py::arg("to"),
             py::arg("weight"),
             py::arg("k"),
             py::arg("kc") = 10,
             R"pbdoc(
            Construct the k-nearest neighbor graph of the indexed MinHash vectors. It will be written to out parameters :obj:`from`, :obj:`to`, and :obj:`weight` as an edge list.
            
            Arguments:
                from (:obj:`VectorUint`): A vector to which the ids for the from vertices are written
                to (:obj:`VectorUint`): A vector to which the ids for the to vertices are written
                weight (:obj:`VectorFloat`): A vector to which the edge weights are written
                k (:obj:`int`): The number of nearest neighbors to be retrieved during the construction of the k-nearest neighbor graph

            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
        )pbdoc")
        .def("get_distance", &LSHForest::GetDistance, R"pbdoc(
            Calculate the Jaccard distance between two MinHash vectors.

            Arguments:
                vec_a (:obj:`VectorUint`): A MinHash vector
                vec_b (:obj:`VectorUint`): A MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc")
        .def("get_weighted_distance", &LSHForest::GetWeightedDistance, R"pbdoc(
            Calculate the weighted Jaccard distance between two MinHash vectors.

            Arguments:
                vec_a (:obj:`VectorUint`): A weighted MinHash vector
                vec_b (:obj:`VectorUint`): A weighted MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc")
        .def("get_all_distances", &LSHForest::GetAllDistances, R"pbdoc(
            Calculate the Jaccard distances of a MinHash vector to all indexed MinHash vectors.

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
            
            Returns:
                :obj:`List` of :obj:`float`: The Jaccard distances
        )pbdoc")
        .def("get_distance_by_id", &LSHForest::GetDistanceById, R"pbdoc(
            Calculate the Jaccard distance between two indexed MinHash vectors.

            Arguments:
                a (:obj:`int`): The id of an indexed MinHash vector
                b (:obj:`int`): The id of an indexed MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc")
        .def("get_weighted_distance_by_id",
             &LSHForest::GetWeightedDistanceById,
             R"pbdoc(
            Calculate the Jaccard distance between two indexed weighted MinHash vectors.

            Arguments:
                a (:obj:`int`): The id of an indexed weighted MinHash vector
                b (:obj:`int`): The id of an indexed weighted MinHash vector

            Returns:
                :obj:`float` The weighted Jaccard distance
        )pbdoc")
        .def("store", &LSHForest::Store, R"pbdoc(
            Serializes the current state of this instance of :obj:`LSHForest` to the disk in binary format. The index is not serialized and has to be rebuilt after deserialization.
        
            Arguments:
                path (:obj:`str`): The path to which to searialize the file
        )pbdoc")
        .def("restore", &LSHForest::Restore, R"pbdoc(
            Deserializes a previously serialized (using :obj:`store()`) state into this instance of :obj:`LSHForest` and recreates the index.
        
            Arguments:
                path (:obj:`str`): The path to the file which is deserialized
        )pbdoc")
        .def("size", &LSHForest::size, R"pbdoc(
            Returns the number of MinHash vectors in this LSHForest instance.

            Returns:
                :obj:`int`: The number of MinHash vectors
        )pbdoc")
        .def("get_hash", &LSHForest::GetHash, R"pbdoc(
            Retrieve the MinHash vector of an indexed entry given its index. The index is defined by order of insertion.

            Arguments:
                a (:obj:`int`): The id of an indexed MinHash vector

            Returns:
                :obj:`VectorUint` The MinHash vector
        )pbdoc")
        .def("clear", &LSHForest::Clear, R"pbdoc(
            Clears all the added data and computed indices from this :obj:`LSHForest` instance.
        )pbdoc");

    py::class_<TestSub>(m, "TestSub", R"pbdoc(
        A test.
    )pbdoc");

    py::class_<PyMinhash>(m, "Minhash", R"pbdoc(
        A generator for MinHash vectors that supports binary, indexed, string and also :obj:`int` and :obj:`float` weighted vectors as input.
    )pbdoc")
        .def(py::init<unsigned int, unsigned int, unsigned int>(),
             py::arg("d") = 128,
             py::arg("seed") = 42,
             py::arg("sample_size") = 128,
             R"pbdoc(
            Constructor for the class :obj:`Minhash`.

            Keyword Arguments:
                d (:obj:`int`): The number of permutations used for hashing
                seed (:obj:`int`): The seed used for the random number generator(s)
                sample_size (:obj:`int`): The sample size when generating a weighted MinHash
        )pbdoc")
        .def("from_binary_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(py::list &)>(&PyMinhash::FromBinaryArray), R"pbdoc(
            Create a MinHash vector from a binary array.

            Arguments:
                vec (:obj:`List`): A Python list containing binary values
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("from_binary_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(std::vector<uint8_t> &)>(&PyMinhash::FromBinaryArray), R"pbdoc(
            Create a MinHash vector from a binary array.

            Arguments:
                vec (:obj:`List`): A Python list containing binary values
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_binary_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::list &)>(&PyMinhash::BatchFromBinaryArray), R"pbdoc(
            Create MinHash vectors from binary arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`List`): A list of lists containing binary values
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_binary_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::array_t<uint8_t> &)>(&PyMinhash::BatchFromBinaryArray), R"pbdoc(
         py::overload_cast<py::array_t<uint8_t>&>(&PyMinhash::BatchFromBinaryArray), R"pbdoc(
            Create MinHash vectors from binary arrays (parallelized).

            Arguments:
                vec (:obj:`Array`): A 2D array containing binary values
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_binary_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(std::vector<std::vector<uint8_t>> &)>(&PyMinhash::BatchFromBinaryArray), R"pbdoc(
            Create MinHash vectors from binary arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`VectorUchar`): A list of vectors containing binary values
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("from_sparse_binary_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(py::list &)>(&PyMinhash::FromSparseBinaryArray), R"pbdoc(
            Create a MinHash vector from a sparse binary array.

            Arguments:
                vec (:obj:`List`): A Python list containing indices of ones in a binary array
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("from_sparse_binary_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(std::vector<uint32_t> &)>(&PyMinhash::FromSparseBinaryArray), R"pbdoc(
            Create a MinHash vector from a sparse binary array.

            Arguments:
                vec (:obj:`VectorUint`): A Python list containing indices of ones in a binary array
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_sparse_binary_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::list &)>(&PyMinhash::BatchFromSparseBinaryArray), R"pbdoc(
         R"pbdoc(
            Create MinHash vectors from sparse binary arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`List`): A list of Python lists containing indices of ones in a binary array
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_sparse_binary_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::array_t<uint32_t> &)>(&PyMinhash::BatchFromSparseBinaryArray), R"pbdoc(
         py::overload_cast<>(&PyMinhash::BatchFromSparseBinaryArray),
         R"pbdoc(
            Create MinHash vectors from sparse binary arrays (parallelized).

            Arguments:
                vec (:obj:`Array`): A 2D array containing indices of ones in a binary array
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_sparse_binary_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(std::vector<std::vector<uint32_t>> &)>(&PyMinhash::BatchFromSparseBinaryArray), R"pbdoc(
         R"pbdoc(
            Create MinHash vectors from sparse binary arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`VectorUint`): A list of vectors containing indices of ones in a binary array
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("from_string_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(std::vector<std::string> &)>(&PyMinhash::FromStringArray), R"pbdoc(
            Create a MinHash vector from a string array.

            Arguments:
                vec (:obj:`List` of :obj:`str`): A vector containing strings
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_string_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(std::vector<std::vector<std::string>> &)>(&PyMinhash::BatchFromStringArray), R"pbdoc(
            Create MinHash vectors from string arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`List` of :obj:`str`): A list of list of strings
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("from_weight_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(py::list &, const std::string &)>(&PyMinhash::FromWeightArray),
             py::arg("vec"),
             py::arg("method") = "ICWS",
             R"pbdoc(
            Create a MinHash vector from a :obj:`List`.

            Arguments:
                vec (:obj:`List`): A Python list containing :obj:`float` values
            
            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("from_weight_array",
             static_cast<std::vector<uint32_t> (PyMinhash::*)(std::vector<float> &, const std::string &)>(&PyMinhash::FromWeightArray),
             py::arg("vec"),
             py::arg("method") = "ICWS",
             R"pbdoc(
            Create a MinHash vector from a :obj:`float` array.

            Arguments:
                vec (:obj:`VectorFloat`): A vector containing :obj:`float` values
            
            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_weight_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::list &, const std::string &)>(&PyMinhash::BatchFromWeightArray),
             py::arg("vecs"),
             py::arg("method") = "ICWS",
             R"pbdoc(
            Create MinHash vectors from :obj:`float` arrays (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`List`): A list of Python lists containing :obj:`float` values

            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_weight_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::array_t<float, 32> &, const std::string &)>(&PyMinhash::BatchFromWeightArray),
             py::arg("vecs"),
             py::arg("method") = "ICWS",
             R"pbdoc(
            Create MinHash vectors from :obj:`float` arrays (parallelized).

            Arguments:
                vecs (:obj:`Array`): A 2D array containing :obj:`float` values

            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_weight_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(py::array_t<double> &, const std::string &)>(&PyMinhash::BatchFromWeightArray),
             py::arg("vecs"),
             py::arg("method") = "ICWS",
             R"pbdoc(
            Create MinHash vectors from :obj:`float` arrays (parallelized).

            Arguments:
                vecs (:obj:`Array`): A 2D array containing :obj:`float` values

            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_weight_array",
             static_cast<std::vector<std::vector<uint32_t>> (PyMinhash::*)(std::vector<std::vector<float>> &, const std::string &)>(&PyMinhash::BatchFromWeightArray),
             py::arg("vecs"),
             py::arg("method") = "ICWS",
             R"pbdoc(
            Create MinHash vectors from :obj:`float` arrays (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`VectorFloat`): A list of vectors containing :obj:`float` values

            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_int_weight_array",
             &PyMinhash::BatchFromIntWeightArray,
             py::arg("vecs"),
             py::arg("divide_by") = 0,
             R"pbdoc(
            Create MinHash vectors from :obj:`int` arrays, where entries are weights rather than indices of ones (parallelized).

            Arguments:
                vecs (:obj:`List` of :obj:`VectorUint`): A list of vectors containing :obj:`int` values
                divide_by (:obj:`int`): A integer by which each value of each vector is divided. Information is lost, but the running time will be lowered. Faster if :obj:`divide_by` is a power of two.
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("get_distance", &PyMinhash::GetDistance, R"pbdoc(
            Calculate the Jaccard distance between two MinHash vectors.

            Arguments:
                vec_a (:obj:`VectorUint`): A MinHash vector
                vec_b (:obj:`VectorUint`): A MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc")
        .def("get_weighted_distance", &PyMinhash::GetWeightedDistance, R"pbdoc(
            Calculate the weighted Jaccard distance between two MinHash vectors.

            Arguments:
                vec_a (:obj:`VectorUint`): A weighted MinHash vector
                vec_b (:obj:`VectorUint`): A weighted MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc");

    m.def("get_clusters",
          &GetClusters,
          py::arg("gp"),
          py::arg("classes"),
          R"pbdoc(
        brief Creates clusters from a minimum spanning tree.
        
        Arguments:
            gp (:obj:`GraphProperties`): A GraphProperties object
            classes (:obj:`VectorUint`): The classes of the vertices
        
        Returns:
          :obj:`float` The Jaccard distance
      )pbdoc");

    m.def("MSDR",
          &MSDR,
          py::arg("gp"),
          R"pbdoc(
        brief Implementation of the MSDR clustering algorithm.
        
        Arguments:
            gp (:obj:`GraphProperties`): A GraphProperties object
        
        Returns:
          :obj:`VectorUint[VectorUint]` The vertex ids divided in clusters
      )pbdoc");
}