/**
 * @file bindings.cc
 * @author Daniel Probst (daenuprobst@gmail.com)
 * @brief Pybind11 bindings for tmap.
 * @version 0.1
 * @date 2019-06-17
 * 
 */

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "layout.hh"
#include "lshforest.hh"
#include "minhash.hh"

using namespace tmap;

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<uint8_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint16_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint32_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint64_t>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);

PYBIND11_MODULE(tmap, m)
{
    py::bind_vector<std::vector<uint8_t>>(m, "VectorUchar", "Unsigned 8-bit int vector.");
    py::bind_vector<std::vector<uint16_t>>(m, "VectorUsmall", "Unsigned 16-bit int vector.");
    py::bind_vector<std::vector<uint32_t>>(m, "VectorUint", "Unsigned 32-bit int vector.");
    py::bind_vector<std::vector<float>>(m, "VectorFloat", "Unsigned 32-bit float vector.");
    py::bind_vector<std::vector<uint64_t>>(m, "VectorUlong", "Unsigned 64-bit int vector.");

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
        .def_readwrite("sl_extra_scaling_steps", &LayoutConfiguration::sl_extra_scaling_steps)
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
        .def_readonly("n_connected_components", &GraphProperties::n_connected_components)
        .def_readonly("n_isolated_vertices", &GraphProperties::n_isolated_vertices)
        .def_readonly("degrees", &GraphProperties::degrees)
        .def_readonly("adjacency_list", &GraphProperties::adjacency_list);

    m.def("layout_from_lsh_forest", &LayoutFromLSHForest, py::arg("lsh_forest"), py::arg("config") = LayoutConfiguration(), py::arg("create_mst") = true, py::arg("clear_lsh_forest") = false, py::arg("weighted") = false, R"pbdoc(
        Create minimum spanning tree or k-nearest neighbor graph coordinates and topology from an :obj:`LSHForest` instance.
        
        Arguments:
            lsh_forest (:obj:`LSHForest`): An :obj:`LSHForest` instance
        
        Keyword Arguments:
            config (:obj:`LayoutConfiguration`, optional): An :obj:`LayoutConfiguration` instance
            create_mst (:obj:`bool`): Whether to create a minimum spanning tree or to return coordinates and topology for the k-nearest neighbor graph
            clear_lsh_forest (:obj:`bool`): Whether to run :obj:`clear()` on the :obj:`LSHForest` instance after k-nearest negihbor graph and MST creation and before layout
            weighted (:obj:`bool`): Whether the MinHash vectors in the :obj:`LSHForest` instance are weighted

        Returns:
            :obj:`Tuple[VectorFloat, VectorFloat, VectorUint, VectorUint, GraphProperties]` The x and y coordinates of the vertices, the ids of the vertices spanning the edges, and information on the graph
    )pbdoc");

    m.def("mst_from_lsh_forest", &MSTFromLSHForest,
          py::arg("lsh_forest"),
          py::arg("k"),
          py::arg("kc") = 10,
          py::arg("weighted") = false, R"pbdoc(
        Create minimum spanning tree topology from an :obj:`LSHForest` instance.
        
        Arguments:
            lsh_forest (:obj:`LSHForest`): An :obj:`LSHForest` instance
            int k (:obj:`int`): The number of nearest neighbors used to create the k-nearest neighbor graph
        
        Keyword Arguments:
            int kc (:obj:`int`): The scalar by which k is multiplied before querying the LSH forest. The results are then ordered decreasing based on linear-scan distances and the top k results returned
            weighted (:obj:`bool`): Whether the MinHash vectors in the :obj:`LSHForest` instance are weighted

        Returns:
            :obj:`Tuple[VectorUint, VectorUint]`: the topology of the minimum spanning tree of the data indexed in the LSH forest
    )pbdoc");

    m.def("layout_from_edge_list", &LayoutFromEdgeList,
          py::arg("vertex_count"), py::arg("edges"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("create_mst") = true, R"pbdoc(
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

    py::class_<LSHForest>(m, "LSHForest", R"pbdoc(
        A LSH forest data structure which incorporates optional linear scan to increase the recovery performance. Most query methods are available in parallelized versions named with a :obj:`batch_` prefix.
    )pbdoc")
        .def(py::init<unsigned int, unsigned int, bool, bool>(), py::arg("d") = 128, py::arg("l") = 8, py::arg("store") = true, py::arg("file_backed") = false, R"pbdoc(
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
        .def("index", &LSHForest::Index, R"pbdoc(
            Index the LSH forest. This has to be run after each time new MinHashes were added.
        )pbdoc")
        .def("is_clean", &LSHForest::IsClean, R"pbdoc(
            Returns a boolean indicating whether or not the LSH forest has been indexed after the last MinHash vector was added.

            Returns:
                :obj:`bool`: :obj:`True` if :obj:`index()` has been run since MinHash vectors have last been added using :obj:`add()` or :obj:`batch_add()`. :obj:`False` otherwise
        )pbdoc")
        .def("query_linear_scan", &LSHForest::QueryLinearScan, py::arg("vec"), py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false, R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query_linear_scan_exclude", &LSHForest::QueryLinearScanExclude, py::arg("vec"), py::arg("k"), py::arg("exclude"), py::arg("kc") = 10, py::arg("weighted") = false, R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                exclude (:obj:`List` of :obj:`VectorUint`) A list of ids of indexed MinHash vectors to be excluded from the search
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query_linear_scan_by_id", &LSHForest::QueryLinearScanById, py::arg("id"), py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false, R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                id (:obj:`int`): The id of an indexed MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("query_linear_scan_exclude_by_id", &LSHForest::QueryLinearScanExcludeById, py::arg("id"), py::arg("k"), py::arg("exclude"), py::arg("kc") = 10, py::arg("weighted") = false, R"pbdoc(
            Query k-nearest neighbors with a LSH forest / linear scan combination. :obj:`k`*:obj:`kc` nearest neighbors are searched for using LSH forest; from these, the :obj:`k` nearest neighbors are retrieved using linear scan.             

            Arguments:
                id (:obj:`int`): The id of an indexed MinHash vector
                k (:obj:`int`): The number of nearest neighbors to be retrieved
            
            Keyword Arguments:
                exclude (:obj:`List` of :obj:`VectorUint`) A list of ids of indexed MinHash vectors to be excluded from the search
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted
            
            Returns:
                :obj:`List` of :obj:`Tuple[float, int]`: The results of the query
        )pbdoc")
        .def("linear_scan", &LSHForest::LinearScan, py::arg("vec"), py::arg("indices"), py::arg("k") = 10, py::arg("weighted") = false, R"pbdoc(
            Query a subset of indexed MinHash vectors using linear scan.

            Arguments:
                vec (:obj:`VectorUint`): The query MinHash vector
                indices (:obj:`VectorUint`) The ids of indexed MinHash vectors that define the subset to be queried

            Keyword arguments:
                k (:obj:`int`): The number of nearest neighbors to be retrieved
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted

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
        .def("get_all_nearest_neighbors", &LSHForest::GetAllNearestNeighbors, py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false, R"pbdoc(
            Get the k-nearest neighbors of all indexed MinHash vectors.

            Arguments:
                k (:obj:`int`): The number of nearest neighbors to be retrieved

            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted

            Returns:
                :obj:`VectorUint` The ids of all k-nearest neighbors
        )pbdoc")
        .def("get_knn_graph", &LSHForest::GetKNNGraph, py::arg("from"), py::arg("to"), py::arg("weight"), py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false, R"pbdoc(
            Construct the k-nearest neighbor graph of the indexed MinHash vectors. It will be written to out parameters :obj:`from`, :obj:`to`, and :obj:`weight` as an edge list.
            
            Arguments:
                from (:obj:`VectorUint`): A vector to which the ids for the from vertices are written
                to (:obj:`VectorUint`): A vector to which the ids for the to vertices are written
                weight (:obj:`VectorFloat`): A vector to which the edge weights are written
                k (:obj:`int`): The number of nearest neighbors to be retrieved during the construction of the k-nearest neighbor graph

            Keyword Arguments:
                kc (:obj:`int`): The factor by which :obj:`k` is multiplied for LSH forest retreival
                weighted (:obj:`bool`): Whether the MinHash vectors in this :obj:`LSHForest` instance are weighted
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
        .def("get_weighted_distance_by_id", &LSHForest::GetWeightedDistanceById, R"pbdoc(
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

    py::class_<Minhash>(m, "Minhash", R"pbdoc(
        A generator for MinHash vectors that supports binary, indexed, string and also :obj:`int` and :obj:`float` weighted vectors as input.
    )pbdoc")
        .def(py::init<unsigned int, unsigned int, unsigned int>(), py::arg("d") = 128, py::arg("seed") = 42, py::arg("sample_size") = 128, R"pbdoc(
            Constructor for the class :obj:`Minhash`.

            Keyword Arguments:
                d (:obj:`int`): The number of permutations used for hashing
                seed (:obj:`int`): The seed used for the random number generator(s)
                sample_size (:obj:`int`): The sample size when generating a weighted MinHash
        )pbdoc")
        .def("from_binary_array", &Minhash::FromBinaryArray, R"pbdoc(
            Create a MinHash vector from a binary array.

            Arguments:
                vec (:obj:`VectorUchar`): A vector containing binary values
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_binary_array", &Minhash::BatchFromBinaryArray, R"pbdoc(
            Create MinHash vectors from binary arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`VectorUchar`): A list of vectors containing binary values
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("from_sparse_binary_array", &Minhash::FromSparseBinaryArray, R"pbdoc(
            Create a MinHash vector from a sparse binary array.

            Arguments:
                vec (:obj:`VectorUint`): A vector containing indices of ones in a binary array
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_sparse_binary_array", &Minhash::BatchFromSparseBinaryArray, R"pbdoc(
            Create MinHash vectors from sparse binary arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`VectorUint`): A list of vectors containing indices of ones in a binary array
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("from_string_array", &Minhash::FromStringArray, R"pbdoc(
            Create a MinHash vector from a string array.

            Arguments:
                vec (:obj:`List` of :obj:`str`): A vector containing strings
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_string_array", &Minhash::BatchFromStringArray, R"pbdoc(
            Create MinHash vectors from string arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`List` of :obj:`str`): A list of list of strings
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("from_weight_array", &Minhash::FromWeightArray, py::arg("vec"), py::arg("method") = "ICWS", R"pbdoc(
            Create a MinHash vector from a :obj:`float` array.

            Arguments:
                vec (:obj:`VectorFloat`): A vector containing :obj:`float` values
            
            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`VectorUint`: A MinHash vector
        )pbdoc")
        .def("batch_from_weight_array", &Minhash::BatchFromWeightArray, py::arg("vecs"), py::arg("method") = "ICWS", R"pbdoc(
            Create MinHash vectors from :obj:`float` arrays (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`VectorFloat`): A list of vectors containing :obj:`float` values

            Keyword Arguments:
                method (:obj:`str`): The weighted hashing method to use (ICWS or I2CWS)
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("batch_from_int_weight_array", &Minhash::BatchFromIntWeightArray, R"pbdoc(
            Create MinHash vectors from :obj:`int` arrays, where entries are weights rather than indices of ones (parallelized).

            Arguments:
                vec (:obj:`List` of :obj:`VectorUint`): A list of vectors containing :obj:`int` values
            
            Returns:
                :obj:`List` of :obj:`VectorUint`: A list of MinHash vectors
        )pbdoc")
        .def("get_distance", &Minhash::GetDistance, R"pbdoc(
            Calculate the Jaccard distance between two MinHash vectors.

            Arguments:
                vec_a (:obj:`VectorUint`): A MinHash vector
                vec_b (:obj:`VectorUint`): A MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc")
        .def("get_weighted_distance", &Minhash::GetWeightedDistance, R"pbdoc(
            Calculate the weighted Jaccard distance between two MinHash vectors.

            Arguments:
                vec_a (:obj:`VectorUint`): A weighted MinHash vector
                vec_b (:obj:`VectorUint`): A weighted MinHash vector

            Returns:
                :obj:`float` The Jaccard distance
        )pbdoc");
}