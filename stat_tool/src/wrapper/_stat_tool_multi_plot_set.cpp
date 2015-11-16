#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_multi_plot_set()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        class ::stat_tool::MultiPlot & (::stat_tool::MultiPlotSet::*method_pointer_4bf7792f9b005726a7e27e69aa348c18)(int) = &::stat_tool::MultiPlotSet::operator[];
        int (::stat_tool::MultiPlotSet::*method_pointer_987c24bb5c8a51ed8fbf5b4ff73f0f9d)() = &::stat_tool::MultiPlotSet::size;
        boost::python::class_< class ::stat_tool::MultiPlotSet, std::shared_ptr< class ::stat_tool::MultiPlotSet > >("MultiPlotSet", boost::python::no_init)
            .def(boost::python::init< int >())
            .def(boost::python::init< int, int >())
            .def("__setitem__", method_pointer_4bf7792f9b005726a7e27e69aa348c18, boost::python::return_internal_reference<>())
            .def("size", method_pointer_987c24bb5c8a51ed8fbf5b4ff73f0f9d)
            .def_readwrite("title", &::stat_tool::MultiPlotSet::title)
            .def_readwrite("border", &::stat_tool::MultiPlotSet::border)
            .def_readwrite("nb_variable", &::stat_tool::MultiPlotSet::nb_variable)
            .def_readwrite("variable_nb_viewpoint", &::stat_tool::MultiPlotSet::variable_nb_viewpoint)
            .def_readwrite("variable", &::stat_tool::MultiPlotSet::variable)
            .def_readwrite("viewpoint", &::stat_tool::MultiPlotSet::viewpoint);
}