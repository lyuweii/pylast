#include "DataFrame.hh"
#include "EventSource.hh"
#include "nanobind/nanobind.h"
#include "nanobind/stl/optional.h"
#include "nanobind/stl/shared_ptr.h"
#include "nanobind/stl/string.h"
#include <arrow/api.h>
#include <arrow/python/pyarrow.h>
#include <arrow/table.h>

namespace nb = nanobind;

PyObject* to_pyarrow(const std::shared_ptr<arrow::Table>& t) {
    // 需要先初始化一次 Arrow <-> Python 桥接
    static bool inited = false;
    if (!inited) {
        arrow::py::import_pyarrow();
        inited = true;
    }
    return arrow::py::wrap_table(t);
}

nb::object to_python(const std::shared_ptr<arrow::Table>& t) {
    return nb::steal(to_pyarrow(t));  // 拿到真正的 pyarrow.Table
}


void bind_dataframe(nb::module_ &m) {
    // 初始化 pyarrow
    //arrow::py::import_pyarrow();

    nb::class_<df::DataTable>(m, "DataTable")
        .def("simulation_table", [](const df::DataTable &dt) {
            return to_python(dt.simulation_table);
        })
        .def("reconstructor_table", [](const df::DataTable &dt) {
            return to_python(dt.reconstructor_table);
        })
        .def("telescope_table",[](const df::DataTable &dt){
            return to_python(dt.telescope_table);
        });

    nb::class_<df::DataFrameMaker>(m, "DataFrame")
        .def(nb::init<EventSource &>(), nb::arg("source"),
             nb::rv_policy::reference_internal)
        .def("__call__", [](df::DataFrameMaker &self) {
            // auto t = self();
            // PyObject *py_table = arrow::py::wrap_table(t);
            // return nb::steal<nb::object>(py_table);
            return self();
        });
}

NB_MODULE(_pylast_dataframe, m) { bind_dataframe(m); }