#include <boost/python.hpp>
#include <stat_tool/stat_tools.h>

void _stat_tool_column_width()
{
        std::string stat_tool_0cdd446515295e8e8373e99f328c3748_name = boost::python::extract< std::string >(boost::python::scope().attr("__name__") + ".stat_tool");
        boost::python::object stat_tool_0cdd446515295e8e8373e99f328c3748_module(boost::python::handle<  >(boost::python::borrowed(PyImport_AddModule(stat_tool_0cdd446515295e8e8373e99f328c3748_name.c_str()))));
        boost::python::scope().attr("stat_tool") = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        boost::python::scope stat_tool_0cdd446515295e8e8373e99f328c3748_scope = stat_tool_0cdd446515295e8e8373e99f328c3748_module;
        int (*function_pointer_98375177f4775501ba8621e3419b6049)(int) = ::stat_tool::column_width;
        boost::python::def("column_width", function_pointer_98375177f4775501ba8621e3419b6049);
        int (*function_pointer_a4ec1929ca0853fbb1225d63f381e9e9)(double, double) = ::stat_tool::column_width;
        boost::python::def("column_width", function_pointer_a4ec1929ca0853fbb1225d63f381e9e9);
        int (*function_pointer_1be9849ad9f95e8cb87cd33a00f772da)(int, int) = ::stat_tool::column_width;
        boost::python::def("column_width", function_pointer_1be9849ad9f95e8cb87cd33a00f772da);
}