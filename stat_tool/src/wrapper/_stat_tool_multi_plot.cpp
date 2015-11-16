#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_multi_plot()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::SinglePlot & (::stat_tool::MultiPlot::*method_pointer_af406258214b5e0d8168d234bb932acb)(int) = &::stat_tool::MultiPlot::operator[];
        void (::stat_tool::MultiPlot::*method_pointer_95aa6653df98514e8e5b7e3f8fbf1b4b)(int) = &::stat_tool::MultiPlot::resize;
        int (::stat_tool::MultiPlot::*method_pointer_eab7c2749bae5da3903dfe580bd77831)() = &::stat_tool::MultiPlot::size;
        boost::python::class_< class ::stat_tool::MultiPlot, std::shared_ptr< class ::stat_tool::MultiPlot > >("MultiPlot", boost::python::no_init)
            .def(boost::python::init< int >())
            .def("__setitem__", method_pointer_af406258214b5e0d8168d234bb932acb, boost::python::return_internal_reference<>())
            .def("resize", method_pointer_95aa6653df98514e8e5b7e3f8fbf1b4b)
            .def("size", method_pointer_eab7c2749bae5da3903dfe580bd77831)
            .def_readwrite("title", &::stat_tool::MultiPlot::title)
            .def_readwrite("xtics", &::stat_tool::MultiPlot::xtics)
            .def_readwrite("ytics", &::stat_tool::MultiPlot::ytics)
            .def_readwrite("xrange", &::stat_tool::MultiPlot::xrange)
            .def_readwrite("yrange", &::stat_tool::MultiPlot::yrange)
            .def_readwrite("xlabel", &::stat_tool::MultiPlot::xlabel)
            .def_readwrite("ylabel", &::stat_tool::MultiPlot::ylabel)
            .def_readwrite("grid", &::stat_tool::MultiPlot::grid)
            .def_readwrite("group", &::stat_tool::MultiPlot::group);
}