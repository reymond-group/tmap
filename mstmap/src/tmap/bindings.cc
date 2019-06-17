#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "layout.hh"
#include "lshforest.hh"
#include "minhash.hh"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<uint8_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint16_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint32_t>);
PYBIND11_MAKE_OPAQUE(std::vector<uint64_t>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);

PYBIND11_MODULE(tmap, m)
{
    py::bind_vector<std::vector<uint8_t>>(m, "VectorUchar");
    py::bind_vector<std::vector<uint16_t>>(m, "VectorUsmall");
    py::bind_vector<std::vector<uint32_t>>(m, "VectorUint");
    py::bind_vector<std::vector<float>>(m, "VectorFloat");
    py::bind_vector<std::vector<uint64_t>>(m, "VectorUlong");

    py::enum_<ScalingType>(m, "ScalingType", py::arithmetic())
        .value("Absolute", ScalingType::Absolute)
        .value("RelativeToAvgLength", ScalingType::RelativeToAvgLength)
        .value("RelativeToDesiredLength", ScalingType::RelativeToDesiredLength)
        .value("RelativeToDrawing", ScalingType::RelativeToDrawing)
        .export_values();

    py::enum_<Placer>(m, "Placer", py::arithmetic())
        .value("Barycenter", Placer::Barycenter)
        .value("Solar", Placer::Solar)
        .value("Circle", Placer::Circle)
        .value("Median", Placer::Median)
        .value("Random", Placer::Random)
        .value("Zero", Placer::Zero)
        .export_values();

    py::enum_<Merger>(m, "Merger", py::arithmetic())
        .value("EdgeCover", Merger::EdgeCover)
        .value("LocalBiconnected", Merger::LocalBiconnected)
        .value("Solar", Merger::Solar)
        .value("IndependentSet", Merger::IndependentSet)
        .export_values();

    py::class_<LayoutConfiguration>(m, "LayoutConfiguration")
        .def(py::init())
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

    py::class_<GraphProperties>(m, "GraphProperties")
        .def(py::init())
        .def_readwrite("mst_weight", &GraphProperties::mst_weight)
        .def_readwrite("n_connected_components", &GraphProperties::n_connected_components)
        .def_readwrite("n_isolated_vertices", &GraphProperties::n_isolated_vertices)
        .def_readwrite("degrees", &GraphProperties::degrees)
        .def_readwrite("adjacency_list", &GraphProperties::adjacency_list);

    m.def("layout_from_lsh_forest", &LayoutFromLSHForest,
          py::arg("lsh_forest"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("create_mst") = true,
          py::arg("clear_lsh_forest") = false,
          py::arg("weighted") = false);

    m.def("mst_from_lsh_forest", &MSTFromLSHForest,
          py::arg("lsh_forest"),
          py::arg("k"),
          py::arg("kc") = 10,
          py::arg("weighted") = false);

    m.def("layout_from_edge_list", &LayoutFromEdgeList,
          py::arg("vertex_count"), py::arg("edges"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("create_mst") = true);

    py::class_<LSHForest>(m, "LSHForest")
        .def(py::init<unsigned int, unsigned int, bool, bool>(), py::arg("d") = 128, py::arg("l") = 8, py::arg("store") = false, py::arg("file_backed") = false)
        .def("add", &LSHForest::Add)
        .def("batch_add", &LSHForest::BatchAdd)
        .def("index", &LSHForest::Index)
        .def("is_clean", &LSHForest::IsClean)
        .def("query_linear_scan", &LSHForest::QueryLinearScan, py::arg("vec"), py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false)
        .def("query_linear_scan_exclude", &LSHForest::QueryLinearScanExclude, py::arg("vec"), py::arg("k"), py::arg("exclude"), py::arg("kc") = 10, py::arg("weighted") = false)
        .def("query_linear_scan_by_id", &LSHForest::QueryLinearScanById, py::arg("id"), py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false)
        .def("query_linear_scan_exclude_by_id", &LSHForest::QueryLinearScanExcludeById, py::arg("id"), py::arg("k"), py::arg("exclude"), py::arg("kc") = 10, py::arg("weighted") = false)
        .def("linear_scan", &LSHForest::LinearScan, py::arg("vec"), py::arg("indices"), py::arg("k") = 0, py::arg("weighted") = false)
        .def("query", &LSHForest::Query)
        .def("queryExclude", &LSHForest::QueryExclude)
        .def("query_by_id", &LSHForest::QueryById)
        .def("query_exclude_by_id", &LSHForest::QueryExcludeById)
        .def("batch_query", &LSHForest::BatchQuery)
        .def("get_all_nearest_neighbors", &LSHForest::GetAllNearestNeighbors, py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false)
        .def("get_knn_graph", &LSHForest::GetKNNGraph, py::arg("from"), py::arg("to"), py::arg("weight"), py::arg("k"), py::arg("kc") = 10, py::arg("weighted") = false)
        .def("get_distance", &LSHForest::GetDistance)
        .def("get_weighted_distance", &LSHForest::GetWeightedDistance)
        .def("get_all_distances", &LSHForest::GetAllDistances)
        .def("get_distance_by_id", &LSHForest::GetDistanceById)
        .def("get_weighted_distance_by_id", &LSHForest::GetWeightedDistanceById)
        .def("store", &LSHForest::Store)
        .def("restore", &LSHForest::Restore)
        .def("size", &LSHForest::size)
        .def("get_hash", &LSHForest::GetHash)
        .def("get_layout", &LSHForest::GetLayout, py::arg("config") = LayoutConfiguration(), py::arg("create_mst") = true, py::arg("mem_dump") = true)
        .def("clear", &LSHForest::Clear);

    py::class_<Minhash>(m, "Minhash")
        .def(py::init<unsigned int, unsigned int, unsigned int>(), py::arg("d") = 128, py::arg("seed") = 42, py::arg("sample_size") = 128)
        .def("from_binary_array", &Minhash::FromBinaryArray)
        .def("batch_from_binary_array", &Minhash::BatchFromBinaryArray)
        .def("from_sparse_binary_array", &Minhash::FromSparseBinaryArray)
        .def("batch_from_sparse_binary_array", &Minhash::BatchFromSparseBinaryArray)
        .def("from_string_array", &Minhash::FromStringArray)
        .def("batch_from_string_array", &Minhash::BatchFromStringArray)
        .def("from_weight_array", &Minhash::FromWeightArray)
        .def("batch_from_weight_array", &Minhash::BatchFromWeightArray)
        .def("batch_from_int_weight_array", &Minhash::BatchFromIntWeightArray)
        .def("get_distance", &Minhash::GetDistance)
        .def("get_weighted_distance", &Minhash::GetWeightedDistance);
}