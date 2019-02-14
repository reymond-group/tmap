#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include "layout.hh"
#include "lshforest.hh"

namespace py = pybind11;

PYBIND11_MAKE_OPAQUE(std::vector<uint32_t>);
PYBIND11_MAKE_OPAQUE(std::vector<float>);

PYBIND11_MODULE(mstmap, m)
{
    py::bind_vector<std::vector<uint32_t>>(m, "VectorUint");
    py::bind_vector<std::vector<float>>(m, "VectorFloat");

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
        .def_readwrite("fme_iterations", &LayoutConfiguration::fme_iterations)
        .def_readwrite("fme_randomize", &LayoutConfiguration::fme_randomize)
        .def_readwrite("fme_threads", &LayoutConfiguration::fme_threads)
        .def_readwrite("sl_repeats", &LayoutConfiguration::sl_repeats)
        .def_readwrite("sl_extra_scaling_steps", &LayoutConfiguration::sl_extra_scaling_steps)
        .def_readwrite("sl_scaling_x", &LayoutConfiguration::sl_scaling_x)
        .def_readwrite("sl_scaling_y", &LayoutConfiguration::sl_scaling_y)
        .def_readwrite("sl_scaling_type", &LayoutConfiguration::sl_scaling_type)
        .def_readwrite("mmm_repeats", &LayoutConfiguration::mmm_repeats)
        .def_readwrite("placer", &LayoutConfiguration::placer)
        .def_readwrite("merger", &LayoutConfiguration::merger)
        .def_readwrite("merger_factor", &LayoutConfiguration::merger_factor)
        .def_readwrite("merger_adjustment", &LayoutConfiguration::merger_adjustment)
        .def("__repr__", &LayoutConfiguration::ToString);

    m.def("layout", &Layout,
          py::arg("vertex_count"), py::arg("from"), py::arg("to"),
          py::arg("config") = LayoutConfiguration(),
          py::arg("weight") = std::vector<float>());

    py::class_<LSHForest>(m, "LSHForest")
        .def(py::init<int, int>(), py::arg("d") = 128, py::arg("l") = 8)
        .def("add", &LSHForest::Add)
        .def("index", &LSHForest::Index)
        .def("is_clean", &LSHForest::IsClean)
        .def("query", &LSHForest::Query);
}